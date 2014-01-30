"""Utility script to clear arrays from scidb after a crash"""
import re

from scidbpy import interface
sdb = interface.SciDBShimInterface('http://localhost:8080')

arrays = sdb.list_arrays()
print len(arrays), "arrays in database"

R = re.compile('py[0-9]{12,14}_[0-9]{5,5}')
auto_arrays = filter(R.match, arrays)
print len(auto_arrays), "auto-generated arrays will be removed"

for arr in auto_arrays:
    print " - removing", arr
    m = sdb.wrap_array(arr)
    m.persistent=False
    del m
