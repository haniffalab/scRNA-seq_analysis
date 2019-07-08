<?php 
	$dr_fname  = "./dr/" . $_GET["dr_name"] . "_coordinates";
	$dr_file = fopen($dr_fname, 'r');
	$dr_array = fgets($dr_file);
	fclose($dr_file);
	echo $dr_array;
?>
