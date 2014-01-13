#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;

require "create_and_save_surface_and_charges.pm"
require "name_generator.pm"
require "read_settings.pm"

my $settings_file_name = $ARGV[0];
my %settings = &ReadSettings($settings_file_name);

# main cycle
for (my $ind = 0; $ind < %settings{"set_size"}; $ind++)
{
	my $current_file_name = &GetMoleculeFullPath(\%settings, $ind);
	my $document = DiscoveryScript::Open( 
	{ Path => $current_file_name } 
	);
	
	my $molecule = $document->Molecules->Item(0);
	&CreateAndSaveSurfaceAndCharges(\$document, \$molecule, %settings, $ind);
}
