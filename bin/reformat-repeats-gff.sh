#!/bin/bash

cat $1 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ =join("\t", @F)."\n"} print $_' > $2
