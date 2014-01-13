#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;
my $i=1;
my $j=1;

for ($i = 1; $i <= 205; $i++) {
	my $a='';
        
        for($j = 1; $j <= 3 - int(log($i)/log(10)); $j++)
        {
            $a = join('', $a, '0');
        }

	my $document = DiscoveryScript::Open(
	    { 
		Path      => join('', 'D:\pirim\pirim', $a, $i, '\pirim', $a, $i, '.mol'),
		FormatType => "mol"
	    }
	);
	my $allAtoms = $document->Atoms;
	my $surface =
	  $document->CreateSolidSurface( $allAtoms, Mdm::surfaceStyleSolvent, False,
	    Mdm::surfaceColorByElectrostaticPotential, 1.4 );
	$document->Save(join('','D:\pirim\srf',$i,'.wrl') ,'wrl' ); 
	$document->Save(join('','D:\pirim\srf',$i,'.mol') ,'mol' );
	my $molecule = $document->Molecules->Item(0);
	$document->CalculateCharges();
	my $j=open(my $fh, '>', join('','D:\pirim\ch',$i,'.txt')) ;
	print $fh  join('',$molecule->AtomCount,"\n");
	foreach my $atom (@$allAtoms)
	{
	    print $fh join(' ',$atom->XYZ->X,$atom->XYZ->Y,$atom->XYZ->Z,$atom->Charge,"\n");
	}
	close $fh;
	$document->Close();
}