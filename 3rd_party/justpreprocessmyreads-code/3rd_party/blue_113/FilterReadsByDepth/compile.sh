#!/bin/bash
mcs FilterReadsByDepth/Program.cs CommonCode/MerStrings.cs -out:FilterReadsByDepth.exe
mono --aot -O=all,-shared FilterReadsByDepth.exe

