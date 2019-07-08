<?php 
	$mot_de_pass = $_GET['mot_de_pass'];
	$password = "b3mZ5GVYeRZEvk8ctsEk";
	if ($mot_de_pass != $password){
		exit('Password not correct');
	}

	$category_name  = $_GET["category"];
	$color_key_file = $category_name . "_color_key.csv";
	$color_key_file = "./categories/" . $color_key_file;
	$color_key_file = fopen($color_key_file, "r");
	$line = fgets($color_key_file );
	$data_string = "";
	while (($line = fgets($color_key_file)) !== false){
		$col_key  = explode(",", $line);
		$cell_name = $col_key[0];
		$cell_col  = $col_key[1];
		$button_html = str_replace('cell_name', str_replace('"', "", $cell_name), "<div style='background-color: cell_color;'><input type = 'checkbox' onchange='toggleCategoryAction()' name = 'cellTypeBtn' id='cell_name' checked /><label for = 'cell_name'> cell_name</label></div>");
		$button_html = str_replace('cell_color', str_replace('"', "", $cell_col), $button_html);
		$data_string = $data_string . $button_html;
	}
	fclose($color_key_file);
	
	$categories_colors_file = "./categories/" . $category_name . "_colors";
	$categories_colors_file = fopen($categories_colors_file, "r");
	$categories_colors = fgets($categories_colors_file);
	fclose($categories_colors_file);
	$data_string = $data_string . "&&" . $categories_colors;
	
	$categories_indices_file = "./categories/" . $category_name . "_indices";
	$categories_indices_file = fopen($categories_indices_file, "r");
	$categories_indices = fgets($categories_indices_file);
	fclose($categories_indices_file);
	$data_string = $data_string . "&&" . $categories_indices;
	
	echo $data_string;
?>
