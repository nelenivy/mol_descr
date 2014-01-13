#!/usr/bin/perl -w
#
use strict;

sub ReadSettings
{
	my ($settings_file_name) = @_;
	open(my $settings_file, '<', $settings_file_name);
	my %settings;

	while(my $line = <$settings_file>)
	{
		chomp($line);
		my ($key, $value) = split(/=/, $line);
		$settings{$key} = $value;
	}
	
	close($settings_file);
	return %settings;
}

1;