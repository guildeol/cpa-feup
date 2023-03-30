#!/bin/python3

from os import system

print("Sequential runs:")
for i in range (24, 32):
    print(f"n = {i + 1}")
    system(f'../build/ex2/sieve_seq {i + 1} {10}')

print("Parallel runs:")
for i in range (24, 32):
    print(f"n = {i + 1}")
    system(f'../build/ex2/sieve_parallel {i + 1} {10}')