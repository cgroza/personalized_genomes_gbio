#!/usr/bin/perl

while(<>)
    {
        ($contig, $size, @rest) = split /\t/;
        print "$contig\t$size\n";
    }
