#!/bin/bash
source activate firedpy
echo '{"irods_host": "data.cyverse.org", "irods_port": 1247, "irods_user_name": "$IPLANT_USER", "irods_zone_name": "iplant"}' | envsubst > $HOME/.irods/irods_environment.json
exec /usr/bin/tini -- ttyd tmux new -A -s ttyd bash
