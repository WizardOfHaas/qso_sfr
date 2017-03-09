#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use MongoDB;
use Parallel::ForkManager;

my $forks = 100;
my $pm = Parallel::ForkManager->new($forks);

my $client = MongoDB->connect('mongodb://localhost');
my $collection = $client->ns('qso_morpho.qso_lines_raw');

open my $cfg_file, "<", $ARGV[0];
chomp(my @lines = <$cfg_file>);

my @batches;
my $i = 0;

while(my $d = pop @lines){
	if($batches[$i % $forks]){
		push(@{$batches[$i % $forks]}, $d);
	}else{
		$batches[$i % $forks] = [$d];
	}

	$i++;
}

foreach my $batch (@batches){
	my $pid = $pm->start and next;

	$client->reconnect;

	foreach my $line (@{$batch}){
		my ($PLATEID, $MJD, $FIBER) = split(",", $line);
		my $url = "http://dr12.sdss3.org/spectrumDetail?mjd=$MJD&fiber=$FIBER&plateid=$PLATEID";
		#print $url."\n"; die;

		my $raw_data = `curl -L "$url"`;
		#print $raw_data; die;

		my @t = split(/<tbody>|<\/tbody>/, $raw_data);

		foreach my $entry (split(/<\/tr>/, $t[-2])){
			my ($name) = $entry =~ /<\!\-\- Line Name \-\-><td class\="name">(.*)<\/td>/;
			my ($wavelength) = $entry =~ /<\!\-\- Wavelength \-\-><td>(.*)<\/td>/;
			my ($z) = $entry =~ /<\!\-\- linez \-\-><td>(.*)<\/td>/;
    		my ($sigma) = $entry =~ /<\!\-\- linesigma \-\-><td>(.*)<\/td>/;
    		my ($area) = $entry =~ /<\!\-\- linearea \-\-><td>(.*)<\/td>/;
    		my ($contlevel) = $entry =~ /<\!\-\- linecontlevel \-\-><td>(.*)<\/td>/;

			if($name){
				#Convert wavelength to v, for refrence purposes.
				my $CIV = 1549;
				my $abs_z = ($wavelength/$CIV) - 1;
				my $RC = (1 + $z)/(1 + $abs_z);
				my $betaC = (($RC**2) - 1)/(($RC**2) + 1);
				my $beta = $betaC * (-30000);

				my $data = {
					name => $name,
					wavelength => 0+ $wavelength,
					z => ($z ? 0+ $z : 0),
					sigma => 0+ $sigma,
					area => 0+ $area,
					contlevel => 0+ $contlevel,
					SID => (0+ $PLATEID)."-".(0+ $MJD)."-".(0+ $FIBER),
					beta => 0+ $beta
				};

				$collection->insert_one($data);

			#print Dumper $data;
			}
		}
	}

	$pm->finish;
}

$pm->wait_all_children;
