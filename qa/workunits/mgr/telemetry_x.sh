#!/bin/sh -ex

# Make sure we're still opted-in

LAST_OPT_REVISION=$(ceph config get mgr mgr/telemetry/last_opt_revision)
DESCRIPTION=$(ceph config get mgr mgr/telemetry/description)
if [ $LAST_OPT_REVISION != 3 ]; then
    echo "Expected last_opt_revision to be 3; got $LAST_OPT_REVISION"
    exit 1
fi

# Check that warning in the health

CEPH_STATUS=$(ceph -s)
if [[ "$CEPH_STATUS" != *"Telemetry requires re-opt-in"* ]]; then
    echo "Expected ceph status to contain re-opt-in warning; got $CEPH_STATUS"
    exit 1
fi

# Re-opt-in

ceph telemetry on --license sharing-1-0

# Wait a bit, then that the warning is gone

sleep 10
CEPH_STATUS=$(ceph -s)
if [[ "$CEPH_STATUS" == *"Telemetry requires re-opt-in"* ]]; then
    echo "Expected ceph status to not contain re-opt-in warning; got $CEPH_STATUS"
    exit 1
fi

# List channels (for debugging)

ceph telemetry channel ls

# Check that the new collections are there

COLLECTIONS=$(ceph telemetry collection ls)
if [[ "$COLLECTIONS" != *"basic_mds_metadata"* ]]; then
    echo "Expected to have the basic_mds_metadata collection; got $COLLECTIONS"
    exit 1
elif [[ "$COLLECTIONS" != *"basic_pool_usage"* ]]; then
    echo "Expected to have the basic_pool_usage collection; got $COLLECTIONS"
    exit 1
elif [[ "$COLLECTIONS" != *"basic_rook_v01"* ]]; then
    echo "Expected to have the basic_rook_v01 collection; got $COLLECTIONS"
    exit 1
elif [[ "$COLLECTIONS" != *"basic_usage_by_class"* ]]; then
    echo "Expected to have the basic_usage_by_class collection; got $COLLECTIONS"
    exit 1
elif [[ "$COLLECTIONS" != *"perf_perf"* ]]; then
    echo "Expected to have the perf_perf collection; got $COLLECTIONS"
    exit 1
fi

# Display reports

ceph telemetry show
ceph telemetry show-device

# Opt out (so we don't keep sending reports)

ceph telemetry off

echo OK
