#!/usr/bin/python3

import random
_b32_symbols = list('0123456789abcdefghjkmnpqrstvwxyz')
def urlb32_encode(i, zeropad=0):
    if not isinstance(i, int):
        raise TypeError('first input must be integer type')
    output = []
    while i > 0:
        digit = i % 32
        i = i // 32
        output.append(_b32_symbols[digit])
    if zeropad > len(output):
        output.extend([ '0' for i in range(zeropad - len(output)) ])
    output2 = []
    for idx in range(len(output)):
        if idx > 0 and idx % 3 == 0:
            output2.append('-')
        output2.append(output[idx])
    output2.reverse()
    return ''.join(output2)

def get_batch_id(): 
    max_int = 32**6 - 1   # 1,073,741,823
    return urlb32_encode(random.randint(0, max_int))


