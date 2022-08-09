#!/bin/sh -ex

# Set up contact info

ceph config set mgr mgr/telemetry/contact "ceph-org"
ceph config set mgr mgr/telemetry/description "upgrade test cluster"
ceph config set mgr mgr/telemetry/channel_ident true

CONTACT=$(ceph config get mgr mgr/telemetry/contact)
DESCRIPTION=$(ceph config get mgr mgr/telemetry/description)
CHANNEL_IDENT=$(ceph config get mgr mgr/telemetry/channel_ident)

if [ "$CONTACT" != "ceph-org" ]; then
    echo "Expected contact to be "ceph-org"; got $CONTACT"
    exit 1
elif [ "$DESCRIPTION" != "upgrade test cluster" ]; then
    echo "Expected description to be "upgrade test cluster"; got $DESCRIPTION"
    exit 1
elif [ "$CHANNEL_IDENT" != "true" ]; then
    echo "Expected channel_ident to be true; got $CHANNEL_IDENT"
    exit 1
fi

# Opt-in to telemetry

ceph telemetry on --license sharing-1-0

LAST_OPT_REVISION=$(ceph config get mgr mgr/telemetry/last_opt_revision)
DESCRIPTION=$(ceph config get mgr mgr/telemetry/description)
if [ $LAST_OPT_REVISION != 3 ]; then
    echo "Expected last_opt_revision to be 3; got $LAST_OPT_REVISION"
    exit 1
fi

# Display reports

ceph telemetry show
ceph telemetry show-device

echo OK
