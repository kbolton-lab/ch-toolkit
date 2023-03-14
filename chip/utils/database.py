import os, sys

import chip.utils.logger as log

import redis
import duckdb

def redis_connect(host, port):
    r = redis.Redis(host=host, port=port, db=0)
    return r

def duckdb_connect_rw(db_file, clobber):
    if os.path.exists(db_file) and clobber == True:
        log.logit(f"Deleting existing duckdb file: {db_file}")
        os.unlink(db_file)
    return duckdb_connect(db_file, read_only=False)

def duckdb_connect_ro(db_file):
    return duckdb_connect(db_file, read_only=True)

def duckdb_connect(db_file, read_only=False):
    if not os.path.exists(db_file) and read_only == True:
        msg = f"Could not find duckdb file on file system: {db_file}"
        log.logit(msg, color="red")
        sys.exit(f"[err] {msg}")

    log.logit(f"Connecting to existing duckdb file: {db_file}")
    connection = duckdb.connect(db_file, read_only=read_only)
    return connection
