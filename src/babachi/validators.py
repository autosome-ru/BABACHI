import re

def convert_frac_to_float(string):
    if re.match(r"^[1-9]+[0-9]*/[1-9]+[0-9]*$", string):
        num, denom = string.split('/')
        if int(denom) <= 0:
            return False
        else:
            value = int(num) / int(denom)
    elif re.match(r"^[1-9]+[0-9]*\.[1-9]+[0-9]*$", string):
        try:
            value = float(string)
        except ValueError:
            return False
    elif re.match(r"^[1-9]+[0-9]*$", string):
        try:
            value = int(string)
        except ValueError:
            return False
    else:
        return False
    if value >= 1:
        return value
    else:
        return False


def check_states(string):
    if not string:
        raise ValueError
    string = string.strip().split(',')
    ret_val = list(map(convert_frac_to_float, string))
    if not all(ret_val):
        raise ValueError
    else:
        return ret_val


def check_samples(string):
    if not string:
        raise ValueError
    string = string.strip().split(',')
    int_conv = []
    for sample in string:
        try:
            int_elem = int(sample)
            if int_elem >= 0:
                int_conv.append(int_elem)
        except ValueError:
            pass
    if len(int_conv) == len(string):
        return int_conv
    else:
        return string