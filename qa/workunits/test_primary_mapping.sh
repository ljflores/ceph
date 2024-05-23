#!/bin/bash

pgid=$(ceph pg dump pgs_brief -f json | jq '.pg_stats[0]' | jq -r '.pgid')
new_prim=$(ceph pg dump pgs_brief -f json | jq '.pg_stats[0]' | jq -r '.up[-1]')

ceph osd set-require-min-compat-client squid
ceph osd pg-upmap-primary "${pgid}" "${new_prim}"
