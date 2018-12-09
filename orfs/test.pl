#!/usr/bin/perl

use strict;
use warnings;

use Tree::Interval::Fast;

my $intervals = Tree::Interval::Fast->new();
my $interval = Tree::Interval::Fast::Interval->new(2, 1, 1);
$intervals->insert($interval);
printf $interval->low;

#$interval = Tree::Interval::Fast::Interval->new(2, 1, 1);
#$intervals->insert($interval);

print $intervals->find(1, 3)->low;
