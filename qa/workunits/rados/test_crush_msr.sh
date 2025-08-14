#!/bin/bash
set -ex

vstart=0
if [ "$1" = "--vstart" ]; then
    vstart=1
fi

ceph="ceph"
if [ $vstart -eq 1 ]; then
    ceph="./bin/ceph"
fi

"$ceph" osd erasure-code-profile set ec84-profile-invalid k=8 m=4 crush-failure-domain=host crush-osds-per-failure-domain=3 crush-num-failure-domains=6

"$ceph" osd set-require-min-compat-client squid

"$ceph" osd pool create ec_test_pool erasure ec84-profile-invalid
