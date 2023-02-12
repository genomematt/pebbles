#!/usr/bin/env bash
minimap2 -x sr -a --MD ref.fa query.fa -o map.sam
samtools view --bam -o map.bam map.sam
#samtools view --cram -O CRAM -o map.cram -T ref.fa map.sam