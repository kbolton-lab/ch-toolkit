import redis

def redis_connect(host, port):
    r = redis.Redis(host=host, port=port, db=0)
    return r
