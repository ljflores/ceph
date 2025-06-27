#!/usr/bin/env bash
set -ex

# Pass in the --vstart flag if you want to run this test locally.
# For example:
#   cd ceph/build
#   ../qa/workunits/rados/test.sh --vstart

parallel=1
[ "$1" = "--serial" ] && parallel=0

# let crimson run in serial mode
crimson=0
[ "$1" = "--crimson" ] && parallel=0 && crimson=1

color=""
[ -t 1 ] && color="--gtest_color=yes"

vstart=0
[ "$1" = "--vstart" ] && vstart=1

function cleanup() {
    pkill -P $$ || true
}
trap cleanup EXIT ERR HUP INT QUIT

GTEST_OUTPUT_DIR=${TESTDIR:-$(mktemp -d)}/archive/unit_test_xml_report
mkdir -p $GTEST_OUTPUT_DIR

declare -A pids

if [ $vstart -eq 1 ]; then
    # compile ceph_test_rados targets
    for f in \
        api_aio api_aio_pp \
	api_io api_io_pp \
	api_asio api_list \
	api_lock api_lock_pp \
	api_misc api_misc_pp \
	api_tier_pp \
	api_pool \
	api_snapshots api_snapshots_pp \
	api_stat api_stat_pp \
	api_watch_notify api_watch_notify_pp \
	api_cmd api_cmd_pp \
	api_service api_service_pp \
	api_c_write_operations \
	api_c_read_operations \
	list_parallel \
	open_pools_parallel \
	delete_pools_parallel
    do
        ninja -j$(nproc) ceph_test_rados_$f
    done

    # compile ceph_test_neorados targets
    for f in \
        cls cmd handler_error io ec_io list ec_list misc pool read_operations snapshots \
        watch_notify write_operations
    do
        ninja -j$(nproc) ceph_test_neorados_$f
    done

    echo "Setting up a test cluster..."
    ninja -j$(nproc) vstart
    ../src/vstart.sh --debug --new -x --localhost --bluestore
fi

# run ceph_test_rados tests
for f in \
    api_aio api_aio_pp \
    api_io api_io_pp \
    api_asio api_list \
    api_lock api_lock_pp \
    api_misc api_misc_pp \
    api_tier_pp \
    api_pool \
    api_snapshots api_snapshots_pp \
    api_stat api_stat_pp \
    api_watch_notify api_watch_notify_pp \
    api_cmd api_cmd_pp \
    api_service api_service_pp \
    api_c_write_operations \
    api_c_read_operations \
    list_parallel \
    open_pools_parallel \
    delete_pools_parallel
do
    executable="ceph_test_rados_$f"
    if [ $vstart -eq 1 ]; then
        executable="./bin/$executable"
    fi
    if [ $parallel -eq 1 ]; then
	r=`printf '%25s' $f`
	ff=`echo $f | awk '{print $1}'`
	bash -o pipefail -exc "$executable --gtest_output=xml:$GTEST_OUTPUT_DIR/$f.xml $color 2>&1 | tee ceph_test_rados_$ff.log | sed \"s/^/$r: /\"" &
	pid=$!
	echo "test $f on pid $pid"
	pids[$f]=$pid
    else
	$executable
    fi
done

# run ceph_test_neorados tests
for f in \
    cls cmd handler_error io ec_io list ec_list misc pool read_operations snapshots \
    watch_notify write_operations
do
    executable="ceph_test_neorados_$f"
    if [ $vstart -eq 1 ]; then
        executable="./bin/$executable"
    fi
    if [ $parallel -eq 1 ]; then
	r=`printf '%25s' $f`
	ff=`echo $f | awk '{print $1}'`
	bash -o pipefail -exc "$executable --gtest_output=xml:$GTEST_OUTPUT_DIR/neorados_$ff.xml $color 2>&1 | tee ceph_test_neorados_$ff.log | sed \"s/^/$r: /\"" &
	pid=$!
	echo "test $f on pid $pid"
	pids[$f]=$pid
    else
	if [ $crimson -eq 1 ]; then
		if [ $f = "ec_io" ] || [ $f = "ec_list" ]; then
			echo "Skipping EC with Crimson"
			continue
		fi
	fi
	$executable
    fi
done

ret=0
if [ $parallel -eq 1 ]; then
for t in "${!pids[@]}"
do
    # Set timeout values
    max_wait=1800  # 30 minutes
    waited=0
    check_interval=10
    pid=${pids[$t]}
    echo "Waiting for test $t (PID $pid)..."
    # Check in a loop with timeout
    # kill -0 checks if the process is running
    # and 2 >/dev/null suppresses error messages if the process is not found
    while kill -0 $pid 2>/dev/null; do
        sleep $check_interval
        waited=$((waited + check_interval))

        if [ $waited -ge $max_wait ]; then
        # Process timed out
        echo "error in $t ($pid) - TIMED OUT after $max_wait seconds"
        kill -9 $pid 2>/dev/null || true
        ret=1
        break
        fi
    done
    # Only wait after process has ended naturally or been killed
    # We only call wait after determining that the process is no longer running
    # So this won't hang indefinitely like https://tracker.ceph.com/issues/70772
    wait $pid 2>/dev/null || {
        echo "Test $t (PID $pid) failed with non-zero exit status"
        ret=1
    }
done
fi

if [ $vstart -eq 1 ]; then
    echo "Shutting down test cluster..."
    ../src/stop.sh
fi

exit $ret
