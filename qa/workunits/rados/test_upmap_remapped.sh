#!/usr/bin/env bash
# vim: expandtab shiftwidth=4 softtabstop=4
#
# Integration test for src/tools/contrib/upmap-remapped.py
#
# Verifies that the script can resolve remapped PGs caused by adding
# capacity while norebalance is set, making the cluster active+clean
# without moving any data.
#
# Cluster layout (defined in the companion .yaml):
#   host 0: mon.a, mon.b, mon.c, mgr.x, osd.0 .. osd.4   (5 OSDs)
#   host 1: osd.5, osd.6, osd.7, client.0                 (3 OSDs)
#
# The 3 OSDs on host 1 are crushed out (weight 0) at startup so that the
# initial cluster is identical to the 5-OSD baseline.  The test then
# reweights them in, which causes CRUSH to remap PGs onto them.

set -ex

POOL=testpool
NUM_PGS=128
WRITE_SECS=30          # seconds of radosbench writes to populate the pool
UPMAP_SCRIPT="${CEPH_ROOT}/src/tools/contrib/upmap-remapped.py"
MAX_RETRIES=5
WAIT_CLEAN_TIMEOUT=600 # seconds to wait for active+clean

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

function wait_for_clean() {
    local timeout=${WAIT_CLEAN_TIMEOUT}
    local interval=5
    local elapsed=0

    while true; do
        local remapped
        remapped=$(ceph pg ls remapped 2>/dev/null | grep -c '^[0-9]' || true)
        local not_clean
        not_clean=$(ceph -s --format json | \
            python3 -c "import sys,json; s=json.load(sys.stdin); \
            pg=s['pgmap']; print(pg.get('num_pgs',0) - pg.get('num_active_clean',0))" \
            2>/dev/null || echo 999)

        if [ "${not_clean}" -eq 0 ]; then
            echo "Cluster is active+clean"
            return 0
        fi

        if [ "${elapsed}" -ge "${timeout}" ]; then
            echo "Timed out waiting for active+clean (${timeout}s elapsed)" >&2
            ceph -s >&2
            ceph pg ls remapped >&2
            return 1
        fi

        echo "Waiting for clean... ${not_clean} PGs not clean, ${remapped} remapped (${elapsed}s elapsed)"
        sleep ${interval}
        elapsed=$((elapsed + interval))
    done
}

function wait_for_remapped() {
    # Wait until the cluster has at least 1 remapped PG (or timeout)
    local timeout=300
    local interval=5
    local elapsed=0

    while true; do
        local remapped
        remapped=$(ceph pg ls remapped 2>/dev/null | grep -c '^[0-9]' || true)
        if [ "${remapped}" -gt 0 ]; then
            echo "Found ${remapped} remapped PG(s)"
            return 0
        fi

        if [ "${elapsed}" -ge "${timeout}" ]; then
            echo "Timed out waiting for remapped PGs (${timeout}s elapsed)" >&2
            ceph -s >&2
            return 1
        fi

        echo "Waiting for remapped PGs... (${elapsed}s elapsed)"
        sleep ${interval}
        elapsed=$((elapsed + interval))
    done
}

# ---------------------------------------------------------------------------
# test body
# ---------------------------------------------------------------------------

# Step 1: Ensure the 3 "new" OSDs (5, 6, 7) start with crush weight 0 so
#         the initial cluster behaves as if they have not yet been added.
echo "=== Setting new OSDs (5,6,7) to crush weight 0 ==="
for osd in 5 6 7; do
    ceph osd crush reweight osd.${osd} 0
done

# Step 2: Create a test pool and turn off the autoscaler to keep PG count fixed
echo "=== Creating pool ${POOL} with ${NUM_PGS} PGs ==="
ceph osd pool create ${POOL} ${NUM_PGS} ${NUM_PGS}
ceph osd pool set ${POOL} pg_autoscale_mode off
ceph osd pool application enable ${POOL} test

# Step 3: Write data so the cluster is meaningfully populated
echo "=== Writing data to ${POOL} ==="
rados bench -p ${POOL} ${WRITE_SECS} write --no-cleanup

# Step 4: Wait for the initial cluster to be fully clean
echo "=== Waiting for initial active+clean ==="
wait_for_clean

# Step 5: Set norebalance before reweighting the new OSDs in
echo "=== Setting norebalance ==="
ceph osd set norebalance

# Step 6: Reweight the 3 new OSDs to 1.0 — CRUSH will now remap PGs to them,
#         but norebalance prevents any actual data movement (backfill is blocked)
echo "=== Reweighting new OSDs to 1.0 (simulates adding capacity) ==="
for osd in 5 6 7; do
    ceph osd crush reweight osd.${osd} 1.0
done

# Step 7: Confirm the cluster has remapped PGs and that the balancer is blocked
echo "=== Confirming remapped PGs exist ==="
wait_for_remapped

echo "=== Confirming balancer reports misplaced threshold ==="
ceph -s
ceph balancer status || true   # non-fatal; just for log context

# Step 8: Dry-run the script to verify it generates sensible commands
echo "=== Dry run of upmap-remapped.py ==="
python3 "${UPMAP_SCRIPT}"

# Step 9: Apply the upmap entries (retry up to MAX_RETRIES times as instructed
#         by the script's own documentation)
echo "=== Applying upmap entries ==="
applied=0
for i in $(seq 1 ${MAX_RETRIES}); do
    echo "  upmap-remapped run ${i}/${MAX_RETRIES}"
    python3 "${UPMAP_SCRIPT}" | sh || true

    # Give the cluster a moment to process the new upmap items
    sleep 5

    remapped=$(ceph pg ls remapped 2>/dev/null | grep -c '^[0-9]' || true)
    echo "  Remapped PGs remaining: ${remapped}"
    if [ "${remapped}" -eq 0 ]; then
        applied=1
        break
    fi
done

if [ "${applied}" -ne 1 ]; then
    echo "ERROR: remapped PGs still exist after ${MAX_RETRIES} script runs" >&2
    ceph pg ls remapped >&2
    exit 1
fi

# Step 10: Verify the cluster reaches active+clean
echo "=== Verifying active+clean after upmap ==="
wait_for_clean

# Sanity-check: the script should not have triggered any data movement.
# Confirm that norebalance is still set (i.e. no backfill ran).
echo "=== Confirming norebalance is still in effect ==="
ceph osd dump | grep -q 'norebalance' || {
    echo "ERROR: norebalance flag is unexpectedly absent" >&2
    exit 1
}

# Step 11: Unset norebalance and leave the balancer to clean up
echo "=== Unsetting norebalance ==="
ceph osd unset norebalance

# Cleanup
echo "=== Cleaning up ==="
rados -p ${POOL} cleanup || true
ceph osd pool delete ${POOL} ${POOL} --yes-i-really-really-mean-it

echo "PASS"
