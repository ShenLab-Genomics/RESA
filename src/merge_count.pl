#!/usr/bin/perl -w
#
$file=$ARGV[0];
open(F,"$file")||die;
while(<F>){
	(@arr)=split /\s+/,$_;
	$var = $arr[1]."_".$arr[2];
	$saw{$var}=1;
	next;
}
print scalar(keys %saw),"\n";
