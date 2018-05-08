#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:33:44 2018

@author: rkoch
"""

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def uniquify(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

