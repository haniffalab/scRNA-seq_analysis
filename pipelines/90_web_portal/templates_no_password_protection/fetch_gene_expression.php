<?php 
	$gene_name  = $_GET["gene_name"];
	$gene_expression_file = './genes/' . $gene_name;
	$gene_expression_file = fopen($gene_expression_file, 'r');
	$gene_expression = fgets($gene_expression_file);
	echo $gene_expression;
?>