#!/bin/bash
mcs MergeCBTFiles/*.cs -out:MergeCBTFiles.exe
mono --aot -O=all,-shared MergeCBTFiles.exe
