import sys, time, datetime

from clint.textui import puts_err, indent, colored

available_colors = {'red', 'green', 'yellow', 'blue', 'black', 'magenta', 'cyan', 'white'}

def colorize(msg, color):
    function = getattr(colored, color)
    formatted_msg = function(msg)
    return formatted_msg

def logit(msg, color=None):
    ts = time.strftime("[ %Y-%m-%d %T ]", datetime.datetime.now().timetuple())
    fullmsg = "{} {}".format(ts, msg)
    formatted_msg = colorize(fullmsg, color) if color in available_colors else fullmsg
    puts_err(formatted_msg)
    sys.stdout.flush()
    sys.stderr.flush()
