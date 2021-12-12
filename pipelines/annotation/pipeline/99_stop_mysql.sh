#!/bin/bash -l

# Define program name
PROGNAME=$(basename $0)

# Load software
module load funannotate/1.7.3_sing

# Define stop mysqldb
stop_mysqldb() { singularity instance stop mysqldb; }

# Set trap to ensure mysqldb is stopped
trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM

# Define error handler
error_exit()
{
    stop_mysqldb
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

# Set some vars
export SINGULARITY_BINDPATH=/bigdata
export SINGULARITYENV_PASACONF=/rhome/jstajich/pasa.CONFIG.template

# Stop Database
stop_mysqldb

