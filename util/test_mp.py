#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################


# https://www.digitalocean.com/community/tutorials/python-multiprocessing-example

from multiprocess import Process, Queue

def f(q):
    q.put('hello world')

if __name__ == '__main__':
    q = Queue()
    p = Process(target=f, args=[q])
    p.start()
    print (q.get())
    p.join()
