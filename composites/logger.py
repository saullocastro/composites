def warn(msg, level=0, silent=False):
    msg = 'WARNING: ' + msg
    if not silent:
        print('\t'*level + msg)
    return msg

def error(msg, level=0, silent=False):
    msg = 'ERROR: ' + msg
    if not silent:
        print('\t'*level + msg)
    return msg

def msg(msg, level=0, silent=False):
    if not silent:
        print('\t'*level + msg)
    return msg
