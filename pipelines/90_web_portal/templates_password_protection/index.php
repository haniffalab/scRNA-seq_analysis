<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Gene expression web portal</title>
  <meta name="description" content="Gene expression web portal for single cell RNA sequencing data made for the Human Cell Atlas at Newcastle University">
  <meta name="author" content="Dorin-Mirel Popescu">
    <style>
body {
  font-family: Avenir, Arial, sans-serif;
}
   #div_left {
    float: left;
}

#canvasdiv {
	overflow-x: scroll;
	overflow-y: scroll;
	max-height: 90vh;
}
#typesControlPanel {
 height: 25em;
 overflow-y: auto;
}

#cursorText{
    position:absolute;
    background-color: black;
    color : white;
}

  </style>
</head>

<body>
	<?php
	$password = 'cucurucu';
	
		echo "<script type = 'text/javascript'>";
		echo "var password = prompt('Introduce password:')";
		echo "</script>";
 	?>

  <div id = "div_left">
  	<form>
		<fieldset>
			<legend><b>Visualisation options</b></legend>
			<label for = 'particleSizeBar'>Particle size: </label>
			<input type='range' name = 'particleSizeBar' min = .3 max = 7 step=0.1 oninput='setParticleSize(value)' value = 2 /><br/>
				
			<label for = 'alphaInput'>Transparency: </label>
			<input type='range' name = 'alphaInput' min = 0 max = 1000 oninput='setAlpha(value)' value = 1000 /><br/>
				
			<label for = 'canvasSizeInput'>Canvas size: </label>
			<input type='range' name = 'canvasSizeInput' min = 200 max = 2000 oninput='setCanvasSize(value)' value = 500 /><br/>
		</fieldset>
	</form>
	<form>
		<fieldset>
			<legend><b>Cell types</b></legend>
			<div>
				<div>
					Select category: <select id='categorySelectMenu' onchange = 'updateCategories(value)'>
						<?php
							$categories_file = fopen("./categories/categories_key", "r");
							$categories_content = fgets($categories_file);
							fclose($categories_file);
							$categories = explode(";", $categories_content);
							foreach($categories as $category){
								$category_fields = explode("->", $category);
								$category_name = $category_fields[0];
								$category_value = $category_fields[1];
								echo "<option value='" . $category_value . "'>" . $category_name . "</option>";
							}
						 ?>
					</select>
				</div>
			</div>
			<label for='toggleRadio'><input type='checkbox' name = 'toggleRadio' id='toggleRadio' onchange='toggleAllTypes()' checked />Show all</label>
			<div id='typesControlPanel'>
				
			</div>
		</fieldset>
	</form>
  </div>
  	
  <div id = "div_right">
  <div>
  	<b>Colour by:</b>
  		<label for='colourType_t'><input type='radio' name='colourType' id='colourType_t' onchange='setColourBy(value)' value='cell_type' checked />Cell type</label>
  		<label for='colourType_g'><input type='radio' name='colourType' id='colourType_g' onchange='setColourBy(value)' value='gene_expression' />Gene expression</label>
  	</div>
  	<table id = "geneSelectorTable">
  		<tr>
  			<td>
  				<label for='geneFamilySelector'>Gene family:</label>
  			</td>
  			<td>
  				<select name='geneFamilySelector' id='geneFamilySelector' onchange='selectFamily(value)'>		
  				</select>
  			</td>
  		</tr>
  		<tr>
  			<td>
  				<label for='geneSymbolSelector'>Gene symbol:</label>
  			</td>
  			<td>
  				<input type = 'text' id='genelist_input' name = 'geneSymbolSelector_datalist' list='geneSymbolSelector_datalist' onchange='getGeneExpression(this)'>
  				<datalist id = 'geneSymbolSelector_datalist'>
  					<select onchange = 'getGeneExpression(this)' id ='geneSymbolSelector'></select>
  				</datalist>
  			</td>
  		</tr>
  	</table>
  	<label for = 'bgColorRadio_white'><input id = 'bgColorRadio_white' name = "bgColorRadio" type = "radio" value = 'white' onchange='setBackground(value)' checked/>White background </label>
  	<label for = 'bgColorRadio_dark'><input id = 'bgColorRadio_dark' name = "bgColorRadio" type = "radio" value = 'dark' onchange='setBackground(value)' />Dark background </label>
  	<br/><span id='expression_scale'></span>
  	<div>
  		Choose coordinates: <select onchange = 'getCoordinates(value)'>
  			<?php
  				$dr_file = fopen("./dr/dr_key", "r");
  				$dr_content = fgets($dr_file);
  				fclose($dr_file);
  				$dr_cats = explode(";", $dr_content);
  				foreach($dr_cats as $dr_cat){
  					echo 1;
  					$dr_fields = explode("->", $dr_cat);
  					$dr_name = $dr_fields[0];
  					echo "<option value='" . $dr_name . "'>" . $dr_name . "</option>";
  				}
  			?>
  		</select>
  	</div>
  	<div id='canvasdiv'><canvas id='canvas' width=500 height=500></canvas></div>
  </div>
  
  <script id='vertex-shader' type='x-shader/x-fragment'>
  	attribute vec4 a_Position;
  	attribute vec3 a_Color;
  	uniform float u_basePointSize;
  	uniform float u_Alpha;
  	uniform int u_PaintFeatureScale;
  	varying vec4 v_Color;
  	void main() {
    	gl_Position = a_Position;
    	gl_PointSize = u_basePointSize;
    	if (u_PaintFeatureScale == 0){
    		v_Color = vec4(a_Color, u_Alpha);
    	}
    	else{
    		float r = 0.0;
    		float g = 0.0;
    		float b = 0.0;
    		r = max(0.0, 2.0 * a_Color.r - 1.0);
    		b = max(0.0, 2.0 * (1.0 - a_Color.r) - 1.0);
    		g = 1.0 - 2.0 * abs(a_Color.r - 0.5);
    		v_Color = vec4(r, g, b, u_Alpha);
    	}
  	}
  </script>
  <script id ='fragment-shader' type='x-shader/x-fragment'>
	precision mediump float;
  	varying vec4 v_Color;
  	void main() {
  		float r = 0.0;
  		vec2 cxy = 2.0 * gl_PointCoord - 1.0;
  		r = dot(cxy, cxy);
  		if (r > 1.0){
  			discard;
  		}
    	gl_FragColor = v_Color;
  	}
  </script>
  
  <?php
  	echo "<script type = 'text/javascript'>";
  	$gene_families_file = fopen("./genes/gene_lists", "r");
  	$gene_families_content = fgets($gene_families_file);
  	fclose($gene_families_file);
  	echo $gene_families_content;
  	
  	$dr_file = fopen("./dr/dr_key", "r");
  	$dr_content = fgets($dr_file);
  	fclose($dr_file);
  	$dr_cats = explode(";", $dr_content);
  	$first_dr = explode('->', $dr_cats[0])[0];
  	echo "first_dr = '" . $first_dr . "'";
  	echo "</script>";
   ?>
  
  <script type = 'text/javascript'>
  	
  	var Matrix4 = function(opt_src) {
  		var i, s, d;
  		if (opt_src && typeof opt_src === 'object' && opt_src.hasOwnProperty('elements')) {
    	s = opt_src.elements;
    	d = new Float32Array(16);
    	for (i = 0; i < 16; ++i) {
      		d[i] = s[i];
    	}
    	this.elements = d;
  		} else {
    	this.elements = new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]);
  		}
	};
	
	Matrix4.prototype.setTranslate = function(x, y, z) {
  		var e = this.elements;
  		e[0] = 1;  e[4] = 0;  e[8]  = 0;  e[12] = x;
  		e[1] = 0;  e[5] = 1;  e[9]  = 0;  e[13] = y;
  		e[2] = 0;  e[6] = 0;  e[10] = 1;  e[14] = z;
  		e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
  		return this;
	};
	
	Matrix4.prototype.setLookAt = function(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ) {
  		var e, fx, fy, fz, rlf, sx, sy, sz, rls, ux, uy, uz;

  		fx = centerX - eyeX;
  		fy = centerY - eyeY;
  		fz = centerZ - eyeZ;

  		// Normalize f.
  		rlf = 1 / Math.sqrt(fx*fx + fy*fy + fz*fz);
  		fx *= rlf;
  		fy *= rlf;
  		fz *= rlf;

  		// Calculate cross product of f and up.
  		sx = fy * upZ - fz * upY;
  		sy = fz * upX - fx * upZ;
  		sz = fx * upY - fy * upX;

  		// Normalize s.
  		rls = 1 / Math.sqrt(sx*sx + sy*sy + sz*sz);
  		sx *= rls;
  		sy *= rls;
  		sz *= rls;

  		// Calculate cross product of s and f.
  		ux = sy * fz - sz * fy;
  		uy = sz * fx - sx * fz;
  		uz = sx * fy - sy * fx;

  		// Set to this.
  		e = this.elements;
  		e[0] = sx;
  		e[1] = ux;
  		e[2] = -fx;
  		e[3] = 0;

  		e[4] = sy;
  		e[5] = uy;
  		e[6] = -fy;
  		e[7] = 0;

  		e[8] = sz;
  		e[9] = uz;
  		e[10] = -fz;
  		e[11] = 0;

  		e[12] = 0;
  		e[13] = 0;
  		e[14] = 0;
  		e[15] = 1;

  		// Translate.
  		return this.translate(-eyeX, -eyeY, -eyeZ);
	};
	
	Matrix4.prototype.translate = function(x, y, z) {
  		var e = this.elements;
  		e[12] += e[0] * x + e[4] * y + e[8]  * z;
  		e[13] += e[1] * x + e[5] * y + e[9]  * z;
  		e[14] += e[2] * x + e[6] * y + e[10] * z;
  		e[15] += e[3] * x + e[7] * y + e[11] * z;
  		return this;
	};
	
	Matrix4.prototype.setPerspective = function(fovy, aspect, near, far) {
  		var e, rd, s, ct;

  		if (near === far || aspect === 0) {
    		throw 'null frustum';
  		}
  		if (near <= 0) {
    		throw 'near <= 0';
  		}
  		if (far <= 0) {
    		throw 'far <= 0';
  		}

  		fovy = Math.PI * fovy / 180 / 2;
  		s = Math.sin(fovy);
  		if (s === 0) {
    		throw 'null frustum';
  		}

  		rd = 1 / (far - near);
  		ct = Math.cos(fovy) / s;

  		e = this.elements;

  		e[0]  = ct / aspect;
  		e[1]  = 0;
  		e[2]  = 0;
  		e[3]  = 0;

  		e[4]  = 0;
  		e[5]  = ct;
  		e[6]  = 0;
  		e[7]  = 0;

  		e[8]  = 0;
  		e[9]  = 0;
  		e[10] = -(far + near) * rd;
  		e[11] = -1;

  		e[12] = 0;
  		e[13] = 0;
  		e[14] = -2 * near * far * rd;
  		e[15] = 0;

  		return this;
	};
  </script>
  
  <script type = 'text/javascript'>
  	function initContext(gl){
  		n = buffer_data_array.length / 5
  		var vertexColourBuffer = gl.createBuffer()
  		gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
  	
  		var FSIZE = buffer_data_array.BYTES_PER_ELEMENT;
  	
  		var a_Position = gl.getAttribLocation(gl.program, 'a_Position')
  		gl.vertexAttribPointer(a_Position, 2, gl.FLOAT, false, FSIZE * 5, 0)
  		gl.enableVertexAttribArray(a_Position)
  	
  		var a_Color = gl.getAttribLocation(gl.program, 'a_Color')
  		gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 5, 2 * FSIZE)
  		gl.enableVertexAttribArray(a_Color)
  	
  		u_basePointSize = gl.getUniformLocation(gl.program, 'u_basePointSize')
  		gl.uniform1f(u_basePointSize, particleSize)
  	
  		u_Alpha = gl.getUniformLocation(gl.program, "u_Alpha")
  		gl.uniform1f(u_Alpha, alphaValue)
  		
  		u_PaintFeatureScale = gl.getUniformLocation(gl.program, 'u_PaintFeatureScale')
  		gl.uniform1i(u_PaintFeatureScale, PaintFeatureScale)
  	
  		gl.clearColor(1, 1, 1, 1);
  		if(bg_color == "dark"){
  			gl.clearColor(0, 0, 0, 1)
  		}
  	
  		gl.disable(gl.DEPTH_TEST)
  		gl.enable(gl.BLEND)
  		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
  	
  		gl.clear(gl.COLOR_BUFFER_BIT);
  		return gl
  	}
  	
  	function getContext(canvasWidget){
  		var names = ['webgl', 'experimental-webgl', 'webkit-3d', 'moz-webgl'];
  		for(var i=0; i<names.length; i++){
  			try{
  				var gl = canvasWidget.getContext(names[i])
  			}catch(e){}
  			if(gl){i=names.length}
  		}
  	
  		var vshader = shadersFromScriptElement(gl, 'vertex-shader', gl.VERTEX_SHADER),
  			fshader = shadersFromScriptElement(gl, 'fragment-shader', gl.FRAGMENT_SHADER)
  			program = gl.createProgram();
  		gl.attachShader(program, vshader)
  		gl.attachShader(program, fshader)
  		gl.linkProgram(program)
  		gl.useProgram(program)
  		gl.program = program
  		return gl
  	}
  	
  	function shadersFromScriptElement(gl, ID, type){
  		shaderScript = document.getElementById(ID)
  		var str = ''
  		var k = shaderScript.firstChild;
  		while(k){
  			if (k.nodeType == 3){
  				str += k.textContent;
  			}
  			k = k.nextSibling
  		}
  		var shader = gl.createShader(type)
  		gl.shaderSource(shader, str)
  		gl.compileShader(shader)
  		return shader
  	}
  	
  	function toggleAllTypes(){
  		for (i=0;i<typesControlPanel.childElementCount;i++){
  			typesControlPanel.children[i].children[0].checked = toggleRadio.checked;
  		}
  		updateBuffer()
  		draw()
  	}
  	
  	function updateBuffer(){
  		var buffer_data = [];
  		// first update indices to be used - for this read the category control panel radio buttons
  		current_indices = []
  		for(i=0;i<typesControlPanel.childElementCount;i++){
  			if(typesControlPanel.children[i].children[0].checked){
  				radio_type = typesControlPanel.children[i].children[0].id
  				current_indices = current_indices.concat(type_indices[radio_type])
  			}
  		}
  		// now just populate the buffer_data
  		if(colour_by == 'gene_expression'){
  			current_indices.forEach(function(index, i){
  				buffer_data.push(dr_coordinates[2 * index])
  				buffer_data.push(dr_coordinates[2 * index + 1])
  				buffer_data.push(gene_expression[index])
  				buffer_data.push(gene_expression[index])
  				buffer_data.push(gene_expression[index])
  			})
  		}else{
  			current_indices.forEach(function(index){
  				buffer_data.push(dr_coordinates[2 * index])
  				buffer_data.push(dr_coordinates[2 * index + 1])
  				buffer_data.push(category_type_colors[3 * index])
  				buffer_data.push(category_type_colors[3 * index + 1])
  				buffer_data.push(category_type_colors[3 * index + 2])	
  			})
  		}
  		buffer_data_array = new Float32Array(buffer_data)
  		n = buffer_data_array.length / 5
  	}
  	
  	function draw(){
  		if(bg_color == "white"){
			gl_context.clearColor(1, 1, 1, 1)
		}else{
			gl_context.clearColor(0, 0, 0, 1)
		}
		gl_context.clear(gl_context.COLOR_BUFFER_BIT);
  		if(n > 0){
  			gl_context.bufferData(gl_context.ARRAY_BUFFER, buffer_data_array, gl_context.STATIC_DRAW)
  			gl_context.drawArrays(gl_context.POINTS, 0, n)
  		}
  	}
  	
  	function setParticleSize(value){
  		particleSize = parseInt(value)
  		gl_context.uniform1f(u_basePointSize, particleSize)
  		updateBuffer()
  		draw()
  	}
  	
  	function setAlpha(value){
  		alphaValue = parseInt(value) / 1000
  		gl_context.uniform1f(u_Alpha, alphaValue)
  		updateBuffer()
  		draw()
  	}
  	
  	function setCanvasSize(value){
  		value = parseInt(value)
  		canvas.width = value
  		canvas.height = value
  		gl_context = getContext(canvas)
		gl_context = initContext(gl_context)
		gl_context.viewport(0, 0, canvas.width, canvas.height)
		updateBuffer()
		draw()
  	}
  	
  	function setBackground(value){
  		bg_color = value;
  		draw()
  	}
  	
  	function toggleCategoryAction(){
  		updateBuffer()
  		draw()
  	}
  	
  	function updateCategories(value){
  		var request;
		try {
			request = new XMLHttpRequest();
		}catch(e){
			try{
				request = new ActiveXObject("Msxml2.XMLHTTP");
			}catch(e){
				try{
					request = new ActiveXObject("Microsoft.XMLHTTP");
				}catch(e){
					alert('Your browser is too old. Update your browser!')
				}
			}
		}
		request.onreadystatechange = function(){
			if (request.readyState == 4){
				response = request.responseText
				if(response == 'Password not correct'){document.body.innerHTML = "<h1>Wrong password. Refresh page an try again</h1>"}
				response = response.split('&&');
				typesControlPanel.innerHTML = response[0]
				toggleRadio.checked = true;
				reponse_colors = response[1].split(',');
				reponse_indices = response[2].split(';');
				category_type_colors = []
				type_indices = [];
				reponse_colors.forEach(function(val){
					category_type_colors.push(parseFloat(val));
				})
				reponse_indices.forEach(function(indices){
					try{
						indices = indices.split('->');
						var indices_name   = indices[0],
							indices_values = indices[1].split(',');
							indices_array = [];
						indices_values.forEach(function(val){
							indices_array.push(parseInt(val))
						})
						type_indices[indices_name] = indices_array;
					}catch(e){}
				})
				toggleCategoryAction()
			}
		}
		queryString = "?category=" + value;
		queryString = queryString + "&mot_de_pass=" + password;
		request.open("GET", "fetch_category.php" + queryString, true)
		request.send(null)
  	}
  	
  	function getGeneExpression(caller){
  		value = caller.value;
  		caller_id = caller.id
  		if (caller.id == 'geneSymbolSelector'){
  			genelist_input.value = value
  		}else{
  			geneSymbolSelector.value = value
  		}
  		gene_expression = []
  		try {
			request = new XMLHttpRequest();
		}catch(e){
			try{
				request = new ActiveXObject("Msxml2.XMLHTTP");
			}catch(e){
				try{
					request = new ActiveXObject("Microsoft.XMLHTTP");
				}catch(e){
					alert('Your browser is too old. Update your browser!')
				}
			}
		}
		request.onreadystatechange = function(){
			if (request.readyState == 4){
				response = request.responseText.split(',');
				gene_expression_scale = parseFloat(response[0]);
				var r = response.shift();
				response.forEach(function(val){
					gene_expression.push(parseFloat(val))
				})
				if(colour_by == 'gene_expression'){
					gl_context = getContext(canvas),
					gl_context = initContext(gl_context);
					toggleCategoryAction()
					if(!isNaN(gene_expression_scale)){
						expression_scale.innerHTML = "<canvas id ='scale_canvas' width = 200 height = 30></canvas><i>" + genelist_input.value + '</i>'
						var scale_canvas = document.getElementById('scale_canvas'),
							scale_context = scale_canvas.getContext('2d');
						scale_gradient = scale_context.createLinearGradient(0, 0, 200, 0);
						scale_gradient.addColorStop(0, 'blue');
						scale_gradient.addColorStop(0.5, 'green');
						scale_gradient.addColorStop(1, 'red');
						scale_context.fillStyle = scale_gradient;
						scale_context.fillRect(0, 20, scale_canvas.width, scale_canvas.height)
						scale_context.fillStyle = 'black'
						scale_context.fillText('0', 10, 10)
						scale_context.fillText(parseInt(10 * gene_expression_scale) / 10, 180, 10)
					}else{expression_scale.innerHTML = ""}
					if (geneFamilySelector.value != 'All'){
						if(gene_families[geneFamilySelector.value].indexOf(genelist_input.value) == -1){
							expression_scale.innerHTML = 'Gene entered is not part of selected gene family!'
						}
					}
					if (gene_list.indexOf(genelist_input.value) == -1){
						expression_scale.innerHTML = 'Gene name mistyped or does not exist!'
					}
					if (genelist_input.value == ''){
						expression_scale.innerHTML = "Choose a gene"
					}
					geneSymbolSelector_datalist.value = genelist_input.value
				}
			}
		}
		value = value.replace('/', "___")
		queryString = "?gene_name=" + value;
		request.open("GET", "fetch_gene_expression.php" + queryString, true)
		request.send(null)
  	}
  	
  	function setColourBy(val){
  		colour_by = val;
  		PaintFeatureScale = 1
  		if (colour_by == "cell_type"){
  		 PaintFeatureScale = 0
  		}
		if(colour_by == "gene_expression"){
			getGeneExpression(genelist_input)
		}else{
			gl_context = getContext(canvas),
			gl_context = initContext(gl_context);
			toggleCategoryAction()
			expression_scale.innerHTML = '';
		}	
  	}
  	
  	function getCoordinates(value){
  		dr_coordinates = []
  		try {
			request = new XMLHttpRequest();
		}catch(e){
			try{
				request = new ActiveXObject("Msxml2.XMLHTTP");
			}catch(e){
				try{
					request = new ActiveXObject("Microsoft.XMLHTTP");
				}catch(e){
					alert('Your browser is too old. Update your browser!')
				}
			}
		}
		request.onreadystatechange = function(){
			if (request.readyState == 4){
				response = request.responseText.split(",")
				response.forEach(function(val){dr_coordinates.push(parseFloat(val))})
				if(canvas_init){
  					updateBuffer()
  					draw()
				}
			}
			
		}
		queryString = "?dr_name=" + value;
		request.open("GET", "fetch_dr_coordinates.php" + queryString, true)
		request.send(null)
  	}
  	
  	function selectFamily(value){
  		geneSymbolSelector_innerHTML = "<select onchange = 'getGeneExpression(this)' id ='geneSymbolSelector'>"
  		if (value == 'All'){
  			gene_list.forEach(function(gene_name, i){
  				geneSymbolSelector_innerHTML = geneSymbolSelector_innerHTML + "<option value = '" + gene_name + "'>" + gene_name + "</option>"
  			})
  		}else{
  			family_genes = gene_families[value]
  			family_genes.forEach(function(gene_name, i){
  				geneSymbolSelector_innerHTML = geneSymbolSelector_innerHTML + "<option value = '" + gene_name + "'>" + gene_name + "</option>"
  			})
  		}
  		geneSymbolSelector_innerHTML = geneSymbolSelector_innerHTML + '</select>'
  		geneSymbolSelector.innerHTML = geneSymbolSelector_innerHTML
  		if (canvas_init){getGeneExpression(genelist_input)}
  	}
  	
  	var category_type_colors = [],
  		type_indices = [],
  		gene_expression = [],
  		dr_coordinates = [],
  		categorySelectMenu = document.getElementById('categorySelectMenu'),
  		genelist_input = document.getElementById('genelist_input'),
  		expression_scale = document.getElementById('expression_scale'),
  		canvas             = document.getElementById('canvas'),
  		typesControlPanel  = document.getElementById('typesControlPanel'),
  		toggleRadio        = document.getElementById('toggleRadio'),
  		geneFamilySelector = document.getElementById('geneFamilySelector'),
  		geneSymbolSelector = document.getElementById('geneSymbolSelector'),
  		particleSize = 5,
  		alphaValue   = 1.0,
  		bg_color     = "white",
  		n            = 0,
  		particleSize = 2,
  		PaintFeatureScale = 0,
  		currentMaxExpression = 0,
  		buffer_data_array = null,
  		colour_by = "cell_types",
  		gene_expression_scale = 0,
  		canvas_init = false;
  		
  	// population gene families
  	geneFamilySelector_innerHTML = "<option value = 'All'>All</option>"
  	for(var gene_family_name in gene_families){
  		geneFamilySelector_innerHTML = geneFamilySelector_innerHTML +  "<option value = '" + gene_family_name + "'>" + gene_family_name + "</option>"
  	}
  	geneFamilySelector.innerHTML = geneFamilySelector_innerHTML
  	
  	selectFamily('All')
  	
  	getCoordinates(first_dr)
  	updateCategories(categorySelectMenu.value)
  	updateBuffer()
  		
  	// create the renderer
  	var gl_context = getContext(canvas),
		gl_context = initContext(gl_context);
	
	draw()
	
	canvas_init = true;
	
	var curTxt=document.createElement('div');
    curTxt.id="cursorText";
    document.body.appendChild(curTxt);
	canvas.addEventListener('mousemove', function(event){
	  var canvasRect = canvas.getBoundingClientRect();
	  var selectedIndex = false;
  	  Xc = 2*(event.clientX - canvasRect.x) / canvas.width - 1;
  	  Yc = -2*(event.clientY - canvasRect.y) / canvas.height + 1;
  	  for(var k=0; k<dr_coordinates.length/2;k++){
	     dd = Math.abs(Xc - dr_coordinates[2*k]) + Math.abs(Yc - dr_coordinates[2*k + 1])
	     if (dd < .003){
	       selectedIndex = k
	       break;
	     }
	   }
	   var selectedCellType = ""
	   for (var xy_cell_type in type_indices){
	      if (type_indices[xy_cell_type].includes(selectedIndex)){
	         selectedCellType = xy_cell_type
	         break
	      }
	   }
	   
	   curTxt.innerHTML = selectedCellType
	   curTxt.style.left = event.pageX + 5 + "px"
	   curTxt.style.top  = event.pageY - curTxt.offsetHeight + "px"
	})
	
	// safari does not support datalist
	// see at https://www.w3schools.com/tags/tryit.asp?filename=tryhtml5_datalist
  
  </script>
  <div style="float: clear;""><hr><span style="font-size:0.8em;">This data portal was created using the web_portal tool (<a href="https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data">github link</a>) developed by Dorin-Mirel Popescu</span><hr></div>
</body>
</html>