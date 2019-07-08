#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  23 glorious 2018

@author: Dorin-Mirel Popescu
"""

import sys
args = sys.argv
output_folder = args[1]
load_from = args[2]
lineage_name = args[3]

insert_lineage = lineage_name.replace("@@", " ")

from os.path import join
save_to = join(output_folder, "index.php")
fetch_php = join(output_folder, "fetch_gene_expression.php")

template = """
<!doctype html>

<html lang='en'>
<head>
  <meta charset='utf-8'>

  <title>Trajectory viewer</title>
  <meta name='description' content='Trajectory analysis of insert_lineage.'>
  <meta name='author' content='Dorin-Mirel Popescu'>
  <style type = 'text/css'>

  	#instructionBtn {
  		margin-bottom: 1em;
  		font-size: 1em;
  		background-color: #00bfff;
    	border: none;
    	color: white;
    	padding: 10px 10px;
    	text-align: center;
    	text-decoration: none;
    	cursor: pointer;
  	}

  	#instructions_div {
  		margin:0; padding:0;
  	}
  	body {
  font-family: Avenir, Arial, sans-serif;
}
  </style>
</head>
<body>
    <div id = "Presentation">
        <input id = 'instructionBtn' type = 'button' value = 'Show description and instructions' onclick = 'showHideInstructions(value)' />
        <div id = 'instructionBoard'>
           <div><h2>Trajectory analysis - insert_lineage </h2></div>
           <div>Development and data analysis: Dorin-Mirel Popescu</div>
           <div>
               <b>Instructions: </b>
               <ul>
                   <li>This portal contains 4 option menus and 2 plotting areas;</li>
                   <li>The options menus are: "Visualisation options", "Colour by:", "Gene expression options", "Cell types:";</li>
                   <li>The plotting areas are: 1) the 3D viewer for difussion map coordinates (lower left) and 2) gene expression levels vs pseudotime (lower right and only visible if "Colour by" is set to "Gene")</li>
                   <li>Click on the 3D viewer area to rotate the object;</li>
                   <li>The visualisation options menu allows the user to set graphics parameters for the 3D viewer</li>
                   <li>The "Colour by:" menu is used to set by what criterion to color the data points in the 3D viewer. These are "Cell type", "Pseudotime" and "Gene";</li>
                   <li>Try these now</li>
                   <li>The "Gene expression options" menu is only functional if the option "Gene" is selected in the menu "Colour by:"</li>
                   <li>In the "Gene expression options" menu the user can select the gene which expression is shown both in the 3D viewer and in the gene expression by pseudotime plot;</li>
                   <li>Gene id can be selected from a drop-down menu or manually typed in the text box right of the drop-down menu;</li>
                   <li>Gene expression can be plotted as it is or min-max normalized. The first approach allows the user the gauge the scale of gene expression and the second approch allows for better visualisation of the changes in gene expression with pseudotime</li>
               </ul>
           </div>
        </div>
    </div>
	<table>
		<tr>
			<td style="text-align:left;">
				<form>
					<fieldset>
						<legend><b>Visualisation options</b></legend>
						<label for = 'particleSizeBar'>Particle size: </label>
						<input type='range' id = "particleSizeBar" name = 'particleSizeBar' min = 10 max = 300 onchange='setParticleSize(value)' value = 150 /><br />
				
						<label for = 'alphaInput'>Transparency: </label>
						<input type='range' id = "alphaInput" name = 'alphaInput' min = 0 max = 1000 onchange='setAlpha(value)' value = 1000 /><br />
				
						<label for = 'canvasSizeInput'>Canvas size: </label>
						<input type='range' id = "canvasSizeInput" name = 'canvasSizeInput' min = 200 max = 2000 onchange='setCanvasSize(value)' value = 500 /><br />
				
						<label for = "zoom">Zoom: </label>
						<input type='range' id = "zoom" name = 'zoom' min = 100 max = 1000 onchange='setZoom(value)' value = 400 /><br />
				
						<label for = 'bgInput'>Dark background: </label>
						<input type='radio' id = "bgInput" name = 'bgkInput' onchange='setBackground(value)' value = 'dark' />
						<label for = 'wgInput'>White background: </label>
						<input type='radio' id = "wgInput" name = 'bgkInput' onchange='setBackground(value)' value = 'white' checked />
						<br />
				
						<label for='sliderX'>Slide X: </label>
						<input type='range' id = "sliderX" name='sliderX' min='-100' max='100' onchange='slideOnX(value)' value='0' />
						<label for='sliderY'>Slide Y: </label>
						<input type='range' id = "sliderY" name='sliderY' min='-100' max='100' onchange='slideOnY(value)' value='0' />
						<br />
					</fieldset>
				</form>
			</td>
			<td style='vertical-align: top' rowspan='2'>
			
				<form>
					<fieldset>
						<legend><b>Colour by:</b></legend>
							<label><input type='radio' name=colourType onchange='setColourByType(value)' value='celltype' checked />Cell type</label><br />
							<label><input type='radio' name=colourType onchange='setColourByType(value)' value='pseudotime' />Pseudotime</label><br />
							<label><input type='radio' name=colourType onchange='setColourByType(value)' value='gene' />Gene</label>
					</fieldset>
				</form>
				<br/>
				<form onsubmit="return false">
					<fieldset>
						<legend><b>Gene expression options</b></legend>
						<label for='geneselector'>Chose gene by ID: </label>
						<select id='geneselector' onchange='newGeneCalled()'>
							gene_options_here
						</select>
                           Or type here <input type = "text" id = "enter_gene_manually" size = "25" onchange = "getEnteredGene()" />
						<br/>
						Gene expression as:<br/>
						<label><input type='radio' name='expressionType' value='snn' onchange='setExpressionType(value)'/>Non-normalized gene expression</label><br/>
						<label><input type='radio' name='expressionType' value='sn' onchange='setExpressionType(value)' checked />Minmax normalized gene expression</label><br/>
					</fieldset>
				</form>
				<br />
				<div>
					<fieldset>
						<legend><b>Cell types:</b></legend>
						<label for='toggleRadio'><input type='checkbox' name = 'toggleRadio' id='toggleRadio' onchange='toggleShowTypes()' checked />Show all:</label>
						<form id = 'ControlPanel'>
							radiocommands
						</form>
					</fieldset>
				</div>
                  <div id = "gene_expression_div" style="width:600px;height:300px;"></div>
			</td>
		</tr> 
		<tr>
			<td style='vertical-align: text-top' >
                  <div id = "expression_scale"></div>
				<canvas id='canvas' width=600 height=600></canvas>
			</td>
		</tr>
	</table>
  <script id='vertex-shader' type='x-shader/x-fragment'>
  	attribute vec4 a_Position;
  	attribute vec3 a_Color;
  	uniform mat4 u_ModelMatrix;
  	uniform mat4 u_ViewMatrix;
  	uniform mat4 u_ProjMatrix;
  	uniform float u_basePointSize;
  	uniform float u_Alpha;
  	varying vec4 v_Color;
  	void main() {
  		vec4 cubePos =  u_ProjMatrix * u_ModelMatrix * u_ViewMatrix * a_Position;
  		float currentWidth = 0.0;
  		currentWidth = 3.0 + (u_basePointSize - 3.0) * (1.0 - cubePos.z / cubePos.w) / 2.0;
    	gl_Position = cubePos;
    	gl_PointSize = currentWidth;
    	v_Color = vec4(a_Color, u_Alpha);
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
  		vec2 D = vec2(0.0, 0.0), centers = vec2(.65, .35);
  		float light = 0.0;
  		light = length(centers - gl_PointCoord);
  		light = .1 + .9 * (pow(50.0, -light));
    	gl_FragColor = v_Color * light + (1.0 - light) * vec4(0.0, 0.0, 0.0, 1.0);
  	}
  </script>
  <script type = "text/javascript" src="plotly.min.js"></script>
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
  <script type='text/javascript'>
     function showHideInstructions(val){
			if (val == 'Show description and instructions'){
				instructionBtn.value = 'Hide description and instructions'
				instructionBoard.style.display = 'block'
			}else{
				instructionBtn.value = 'Show description and instructions'
				instructionBoard.style.display = 'none'
			}
		}
    showHideInstructions("Hide description and instructions")
  	function slideOnX(value){
  		Xshift = parseInt(value);
  		modelMatrix.setTranslate(Xshift, Yshift, 0);
  		gl_context.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function slideOnY(value){
  		Yshift = parseInt(value)
  		modelMatrix.setTranslate(Xshift, Yshift, 0);
  		gl_context.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function setColourByType(value){
  		colourKey = value;
  		colourByType()
  	}
  	
  	function colourByType(){
  		if(colourKey == 'celltype'){
              expression_scale.innerHTML = ""
  			colourByCellType()
  		}else if(colourKey == 'pseudotime'){
              expression_scale.innerHTML = "<canvas id ='scale_canvas' width = 200 height = 30></canvas>"
  			colourByPseudotime()
  		}else{
              expression_scale.innerHTML = "<canvas id ='scale_canvas' width = 200 height = 30></canvas>"
  			colourByGene()
  		}
  	}
          

function newGeneCalled(){
  if (colourKey == "gene"){
    expression_scale.innerHTML = "<canvas id ='scale_canvas' width = 200 height = 30></canvas>"
    colourByGene()
  }else{
    gene_expression_div.innerHTML = ""
    current_gene = geneselector.value;
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
        gene_raw_expression = []
        response.forEach(function(val){
          gene_raw_expression.push(parseFloat(val))
        })
        
        if(expressionType == 'snn'){
          var vector = adaptiveMovingAverage(gene_raw_expression)
          var gene_colors = valuesToColours(vector, 0, 6)
          var max_shown_value = maxRawExpression
          plotGeneExpression(vector, max_shown_value, "Non-normalized gene expression")
        }else{
          var vector = adaptiveMovingAverage(gene_raw_expression)
          vector = minMaxNormalization(vector)
          var gene_colors = valuesToColours(vector, 0, 1)
          var max_shown_value = 1
          plotGeneExpression(vector, max_shown_value, "Normalized gene expression")
        }
      }    
    }
    queryString = "?gene_name=" + current_gene;
    request.open("GET", "fetch_gene_expression.php" + queryString, true)
    request.send(null)
  }
}
  	
  	function colourByCellType(){
  		loadBuffer(selectData(), data_buffer)
  		drawBuffers()
  	}
  	
  	function colourByPseudotime(){
         var scale_canvas = document.getElementById('scale_canvas'),
         scale_context = scale_canvas.getContext('2d');
         scale_context.fillStyle = "white"
         scale_context.fillRect(0, 0, scale_canvas.width, scale_canvas.height)
         scale_gradient = scale_context.createLinearGradient(0, 0, 200, 0)
         scale_gradient.addColorStop(0, "blue")
         scale_gradient.addColorStop(.5, "green")
         scale_gradient.addColorStop(1, "red")
         scale_context.fillStyle = scale_gradient
         scale_context.fillRect(0, 20, scale_canvas.width, scale_canvas.height)
         scale_context.textAlign = "left"
         scale_context.fillStyle = "black"
         scale_context.fillText("0.0", 1, 10)
         scale_context.textAlign = "right"
         scale_context.fillText("1.0", 199, 10)
  		current_pseudotime_buffer = new Float32Array(data_buffer.length)
  		current_pseudotime_buffer.set(data_buffer)
  		points_n = data_buffer.length / 6
  		for (i=0;i<points_n;i++){
  			current_pseudotime_buffer[6 * i + 3] = pseudotime_buffer[3*i]
  			current_pseudotime_buffer[6 * i + 4] = pseudotime_buffer[3*i + 1]
  			current_pseudotime_buffer[6 * i + 5] = pseudotime_buffer[3*i + 2]
  		}
  		loadBuffer(selectData(), current_pseudotime_buffer)
  		drawBuffers()
  	}
  	
  	function setExpressionType(value){
  		expressionType = value
  		newGeneCalled()
  	}
      
    function plotGeneExpression(vector, high_range, plot_title){
       var trace1 = {
          x: ordered_pdt,
          y: vector,
          type: "scatter",
          mode: "lines",
          name: plot_title,
          fill: 'tozeroy',
          line: {
             color: 'rgb(219, 64, 82)',
             width: 0.01
          }
       };
       lab_centers_y = []
       for(var k = 0; k < lab_centers.length; k++){
         lab_centers_y.push(high_range)
       }
       var trace2 = {
          x: lab_centers,
          y: lab_centers_y,
          mode: 'text',
          text: cell_labels_unique,
          textposition: 'top',
          type: 'scatter',
          textfont: {
            family: 'sans serif',
            size: 12,
            color: '#1f77b4'
          }
       }
       var layout = {
          width: 600,
          height: 300,
          yaxis: {range: [0, 1.2 * high_range], title: plot_title},
          xaxis: {title: "Pseudo-time"},
          title: geneselector.value,
          showlegend: false
       };
       dat = [trace1, trace2]
       layout["shapes"] = []
       for (var k = 0; k <= lab_transitions.length; k++){
         if (k == 0){
           x0 = 0; x1 = lab_transitions[k];
         }else if (k == lab_transitions.length){
            x0 = lab_transitions[k-1]; x1 = 1.0
         }else {
            x0 = lab_transitions[k-1]; x1 = lab_transitions[k]
         }
         
         trace = {
           "type": "rect",
           "x0": x0,
           "y0": 0,
           "x1": x1,
           "y1": 1.25*high_range,
            'fillcolor': unique_colors[k],
           'line': {
                'color': 'rgba(128, 0, 128, 1)',
                'width': 0
            }
         }
         layout["shapes"].push(trace)
       }
       Plotly.newPlot("gene_expression_div", dat, layout)
    }
  	
  	function colourByGene(){
         gene_expression_div.innerHTML = ""
  		current_gene = geneselector.value;
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
                 gene_raw_expression = []
                 response.forEach(function(val){
    		         gene_raw_expression.push(parseFloat(val))
    		        })
                 
                 if(expressionType == 'snn'){
          			var vector = adaptiveMovingAverage(gene_raw_expression)
          			var gene_colors = valuesToColours(vector, 0, 6)
                       var max_shown_value = maxRawExpression
                       plotGeneExpression(vector, max_shown_value, "Non-normalized gene expression")
          		}else{
          			var vector = adaptiveMovingAverage(gene_raw_expression)
          			vector = minMaxNormalization(vector)
          			var gene_colors = valuesToColours(vector, 0, 1)
                       var max_shown_value = 1
                       plotGeneExpression(vector, max_shown_value, "Normalized gene expression")
          		}
                  var scale_canvas = document.getElementById('scale_canvas'),
                      scale_context = scale_canvas.getContext('2d');
                  scale_context.fillStyle = "white"
                  scale_context.fillRect(0, 0, scale_canvas.width, scale_canvas.height)
                  scale_gradient = scale_context.createLinearGradient(0, 0, 200, 0)
                  scale_gradient.addColorStop(0, "blue")
                  scale_gradient.addColorStop(.5, "green")
                  scale_gradient.addColorStop(1, "red")
                  scale_context.fillStyle = scale_gradient
                  scale_context.fillRect(0, 20, scale_canvas.width, scale_canvas.height)
                  scale_context.fillStyle = "black"
                  scale_context.textAlign = "left"
                  scale_context.fillText("0", 1, 10)
                  scale_context.textAlign = "right"
                  scale_context.fillText(parseInt(10 * max_shown_value) / 10, 199, 10)
                 
          		genecolor_buffer = new Float32Array(data_buffer.length)
          		genecolor_buffer.set(data_buffer)
          		points_n = data_buffer.length / 6
          		for (i=0;i<points_n;i++){
          			genecolor_buffer[6 * i + 3] = gene_colors[3*i]
          			genecolor_buffer[6 * i + 4] = gene_colors[3*i + 1]
          			genecolor_buffer[6 * i + 5] = gene_colors[3*i + 2]
          		}
          		loadBuffer(selectData(), genecolor_buffer)
          		drawBuffers()
             }    
         }
		queryString = "?gene_name=" + current_gene;
		request.open("GET", "fetch_gene_expression.php" + queryString, true)
		request.send(null)
  	}
  	
  	function valuesToColours(vector, minimum, maximum){
  		colours = []
  		range = maximum - minimum;
  		middle = (maximum + minimum) / 2.0;
  		vector.forEach(function(val, i){
  			r = Math.max(0, 2 * (val - minimum) / range - 1)
  			b = Math.max(0, 2 * (maximum - val) / range - 1)
  			g = 1.0 - 2 * Math.abs(val - middle) / range
  			colours = colours.concat([r, g, b])
  		})
  		colours = new Float32Array(colours);
  		return colours;
  	}
          
    
function update3DPlot(){
  if(colourKey == 'celltype'){
    colourByCellType()
  }else if(colourKey == 'pseudotime'){
    colourByPseudotime()
  }else{
    colourByGene3DPlot()
  }
}

function colourByGene3DPlot(){
  if(expressionType == 'snn'){
    var vector = adaptiveMovingAverage(gene_raw_expression)
    var gene_colors = valuesToColours(vector, 0, 6)
    var max_shown_value = maxRawExpression
  }else{
    var vector = adaptiveMovingAverage(gene_raw_expression)
    vector = minMaxNormalization(vector)
    var gene_colors = valuesToColours(vector, 0, 1)
    var max_shown_value = 1
  }
  
  genecolor_buffer = new Float32Array(data_buffer.length)
  genecolor_buffer.set(data_buffer)
  points_n = data_buffer.length / 6
  for (i=0;i<points_n;i++){
    genecolor_buffer[6 * i + 3] = gene_colors[3*i]
    genecolor_buffer[6 * i + 4] = gene_colors[3*i + 1]
    genecolor_buffer[6 * i + 5] = gene_colors[3*i + 2]
  }
  loadBuffer(selectData(), genecolor_buffer)
  drawBuffers()
}
  	
  	function adaptiveMovingAverage(vector){
  		var colours = [],
  			kernel = 10,
  			minim_kernel = 10,
  			range_factor = 5,
  			window = 2 * kernel;
  		for(i=0;i<vector.length;i++){
  			var start_index = Math.max(1, i - kernel),
  				stop_index  = Math.min(vector.length, i + kernel),
  				local_sd = vector.slice(start_index, stop_index);
  			local_mean = local_sd.reduce(function(sum, val){return sum + val}, 0) / local_sd.length;
  			sqDiffs = local_sd.map(function(value){var diff = value - local_mean; return diff*diff});
  			local_sd = Math.sqrt(sqDiffs.reduce(function(sum, val){return sum + val}, 0))
  			local_kernel = minim_kernel + Math.round(range_factor / (local_sd + .1))
  			start_index = Math.max(1, i - local_kernel)
    		stop_index = Math.min(vector.length, i + local_kernel)
    		local_v = vector.slice(start_index, stop_index);
  			colours.push(local_v.reduce(function(sum, val){return sum + val}, 0) / local_v.length)
  		}
  		
  		return colours
  	}
  	
  	function minMaxNormalization(vector){
  		var minim = vector.reduce(function(a, b){return(Math.min(a, b))})
  		vector = vector.map(function(value){return value - minim})
  		var maxim = vector.reduce(function(a, b){return(Math.max(a, b))});
  		vector = vector.map(function(value){return value / maxim})
  		return vector
  	}
  
  	function selectData(){
  		controlPanel = document.getElementById('ControlPanel')
  		controlRadios = controlPanel.elements
  		values = []
  		for(i=0;i<controlRadios.length;i++){
  			if(controlRadios[i].checked){
  				values = values.concat(index_table[controlRadios[i].id])
  			}
  		}
  		new_indices = []
  		for (i=0;i<values.length;i++){
  			v = values[i]
  			new_indices.push(6*v)
  			new_indices.push(6*v+1)
  			new_indices.push(6*v+2)
  			new_indices.push(6*v+3)
  			new_indices.push(6*v+4)
  			new_indices.push(6*v+5)
  		}
  		return new_indices
  	}
  	
  	function loadBuffer(new_indices, data_buffer_from){
  		current_data_buffer = []
  		new_indices.forEach(function(val, i){current_data_buffer.push(data_buffer_from[val])})
  		current_data_buffer = new Float32Array(current_data_buffer)
  		
  		gl_context.bufferData(gl_context.ARRAY_BUFFER, current_data_buffer, gl_context.STATIC_DRAW); // load data to buffer
  		n = current_data_buffer.length / 6
  	}
  	
  	function drawBuffers(){
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function toggleShowTypes(){
  		toggleRadio = document.getElementById('toggleRadio')
  		controlPanel = document.getElementById('ControlPanel')
  		controlRadios = controlPanel.elements
  		for(i=0;i<controlRadios.length;i++){
  			controlRadios[i].checked = toggleRadio.checked
  		}
  		update3DPlot()
  	}
  
  	function setParticleSize(value){
  		particleSize = parseInt(value)
  		gl_context.uniform1f(u_basePointSize, particleSize)
  		colourByType()
  	}
  	
  	function setAlpha(value){
  		alphaValue = parseInt(value) / 1000
  		gl_context.uniform1f(u_Alpha, alphaValue)
  		colourByType()
  	}
  	
  	function setCanvasSize(value){
  		value = parseInt(value)
  		canvas.width = value
  		canvas.height = value
  		gl_context = getContext(canvas)
		gl_context = initContext(gl_context)
		gl_context.viewport(0, 0, canvas.width, canvas.height)
		if(bg_color == "white"){
			gl_context.clearColor(1, 1, 1, 1)
		}else{
			gl_context.clearColor(0, 0, 0, 1)
		}
		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
  	}
  	
  	function setZoom(value){
  		eyeVN = parseInt(value)
  		farField = eyeVN + 100;
  		rotateData(0, 0)
  	}
  	
  	function setBackground(value){
  		if(value == "dark"){
  			gl_context.clearColor(0, 0, 0, 1)
  			bg_color = "dark"
  		}else{
  			gl_context.clearColor(1, 1, 1, 1)
  			bg_color = "white"
  		}
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT)
  		gl_context.drawArrays(gl_context.POINTS, 0, n)
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
  	
  	function initContext(gl){
  		n = current_data_buffer.length / 6
  		var vertexColourBuffer = gl.createBuffer()
  		gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
  		gl.bufferData(gl.ARRAY_BUFFER, current_data_buffer, gl.STATIC_DRAW)
  	
  		var FSIZE = data_buffer.BYTES_PER_ELEMENT;
  	
  		var a_Position = gl.getAttribLocation(gl.program, 'a_Position')
  		gl.vertexAttribPointer(a_Position, 3, gl.FLOAT, false, FSIZE * 6, 0)
  		gl.enableVertexAttribArray(a_Position)
  	
  		var a_Color = gl.getAttribLocation(gl.program, 'a_Color')
  		gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 6, 3 * FSIZE)
  		gl.enableVertexAttribArray(a_Color)
  	
  		u_basePointSize = gl.getUniformLocation(gl.program, 'u_basePointSize')
  		gl.uniform1f(u_basePointSize, particleSize)
  	
  		u_Alpha = gl.getUniformLocation(gl.program, "u_Alpha")
  		gl.uniform1f(u_Alpha, alphaValue)
  		
  		u_ModelMatrix  = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
  		u_ViewMatrix   = gl.getUniformLocation(gl.program, 'u_ViewMatrix');
  		u_ProjMatrix   = gl.getUniformLocation(gl.program, 'u_ProjMatrix');
  		
  		modelMatrix = new Matrix4(); // The model matrix
  		viewMatrix  = new Matrix4();  // The view matrix
  		projMatrix  = new Matrix4();  // The projection matrix
  		
  		modelMatrix.setTranslate(0, 0, 0);  // 
  		viewMatrix.setLookAt(eyeX, eyeY, eyeZ, 0, 0, 0, upX, upY, upZ); // eyeX, eyeY, eyeZ, camX, camY, camZ, upX, upY, upY
  		projMatrix.setPerspective(30, canvas.width/canvas.height, 100, farField); // fov, ratio, near, far
  		// Pass the model, view, and projection matrix to the uniform variable respectively
  		gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  		gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
  		gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
  		
  		gl.clearColor(1, 1, 1, 1); // add ternary conditional
  	
  		gl.enable(gl.DEPTH_TEST)
  		gl.enable(gl.BLEND)
  		gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
  		//gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA)
  	
  		gl.clear(gl.COLOR_BUFFER_BIT);
  		return gl
  	}
      
     lab_centers_here
     
     cell_labels_unique_here
     
     lab_transitions_here
     
     unique_colors_here
  	
  	var canvas = document.getElementById('canvas'),
  		particleSize = 150,
  		alphaValue = 1.0,
  		bg_color = "white",
  		eyeX    = 0.0,
  		eyeY    = 0.0,
  		eyeZ    = 400.0,
  		upX     = 0.0,
  		upY     = 1.0,
  		upZ     = 0.0,
  		eyeVN   = 400.0,
  		farField = 500.0,
  		previousX = null,
  		previousY = null,
  		currentX  = null,
  		currentY  = null,
  		Xshift    = 0,
  		Yshift    = 0,
  		colourKey = 'celltype',
  		expressionType = 'sn',
  		geneselector = document.getElementById('geneselector'),
         enter_gene_manually = document.getElementById('enter_gene_manually');
         all_gene_names = []
    for(var i = 0; i < geneselector.children.length; i++){
      all_gene_names.push(geneselector.children[i].value)
    }
    function getEnteredGene(){
       entered_gene = enter_gene_manually.value
       if (all_gene_names.includes(entered_gene)){
          enter_gene_manually.value = ""
          geneselector.value = entered_gene
          colourByGene()
       }else{
          enter_gene_manually.value = entered_gene + "not found or has not pass the test statistics"
       }
    }
  		
  	data_buffer = new Float32Array([
  		datahere
  	])
  	
  	pseudotime_buffer = new Float32Array([
  		pseudotime_here
  	])
     ordered_pdt = pseudotime_buffer
  	pseudotime_buffer = valuesToColours(pseudotime_buffer, 0.0, 1.0)
  	
  	gene_raw_expression = []
  	current_gene_here
  	var maxRawExpression = maxRawExpression_here
  	
  	index_table = []
  	indiceshere
  	
  	current_data_buffer = data_buffer
  	
  	gl_context = getContext(canvas)
	gl_context = initContext(gl_context)
  	gl_context.drawArrays(gl_context.POINTS, 0, n)
  	
  	function negCrossProduct(vecA, vecB){
  		crossproduct = [ - vecA[1] * vecB[2] + vecA[2] * vecB[1],
  					     - vecA[2] * vecB[0] + vecA[0] * vecB[2],
  					     - vecA[0] * vecB[1] + vecA[1] * vecB[0]
  		]
  		return(crossproduct)
  	}
  	
  	function vectNorm(vector){
  		return(Math.sqrt((vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2])))
  	}
  	
  	function rotateData(hAngle, vAngle){
  		// change vector for very small angles is approximately the cross product of the eye vector and up vector
  		change = negCrossProduct([eyeX, eyeY, eyeZ], [upX, upY, upZ])
  		// normalize the change vector
  		normChange = vectNorm(change)
  		// scale the change vector by the horizontal angle
  		change = [hAngle * change[0]/normChange, hAngle * change[1]/normChange, hAngle * change[2]/normChange]
  		// update the eye vector by adding the change vector
  		eyeX = eyeX - change[0]
  		eyeY = eyeY - change[1]
  		eyeZ = eyeZ - change[2]
  		// renormalize the eye vector, other it will increase with each change (due to approx error)
  		normEye = vectNorm([eyeX, eyeY, eyeZ])
  		eyeX = eyeVN * eyeX /  normEye
  		eyeY = eyeVN * eyeY /  normEye
  		eyeZ = eyeVN * eyeZ /  normEye
  		
  		// get the (eye, up) plane normal
  		planeInvNormal = negCrossProduct([eyeX, eyeY, eyeZ], [upX, upY, upZ])
  		// in the case of vertical angle, the up vector is already the change vector
  		normChange = vectNorm([upX, upY, upZ])
  		change = [vAngle * upX / normChange, vAngle * upY / normChange, vAngle * upZ / normChange]
  		// update the eye vector by adding the change vector
  		eyeX = eyeX + change[0]
  		eyeY = eyeY + change[1]
  		eyeZ = eyeZ + change[2]
  		// renormalize the eye vector, other it will increase with each change (due to approx error)
  		normEye = Math.sqrt((eyeX * eyeX)+(eyeY * eyeY)+(eyeZ * eyeZ))
  		eyeX = eyeVN * eyeX /  normEye
  		eyeY = eyeVN * eyeY /  normEye
  		eyeZ = eyeVN * eyeZ /  normEye
  		// but the up vector needs changing as well
  		newUp = negCrossProduct([eyeX, eyeY, eyeZ], planeInvNormal)
  		newUpNormal = vectNorm(newUp)
  		upX = -newUp[0] / newUpNormal
  		upY = -newUp[1] / newUpNormal
  		upZ = -newUp[2] / newUpNormal
  		
  		gl_context.clear(gl_context.COLOR_BUFFER_BIT);
  		viewMatrix.setLookAt(eyeX, eyeY, eyeZ, 0, 0, 0, upX, upY, upZ);
  		projMatrix.setPerspective(30, canvas.width/canvas.height, 100, farField); 
  		gl_context.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
  		gl_context.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
  		gl_context.drawArrays(gl_context.POINTS, 0, n);
  	}
  	
  	function startRotating(ev){
  		previousX = ev.clientX
  		previousY = ev.clientY
  		canvas.addEventListener('mousemove', rotateEvent)
  		canvas.addEventListener('mouseup', stopRotation)
  		canvas.addEventListener('mouseout', stopRotation)
  	}
  	
  	function stopRotation(ev){
  		canvas.removeEventListener('mousemove', rotateEvent)
  		canvas.removeEventListener('mouseup', stopRotation)
  		canvas.removeEventListener('mouseout', stopRotation)
  	}
  	
  	function rotateEvent(ev){
  		currentX = ev.clientX
  		currentY = ev.clientY
  		var dX = currentX - previousX,
  			dY = currentY - previousY;
  		rotateData(2.0 * dX, 2.0 * dY)
  		previousX = currentX;
  		previousY = currentY;
  	}
      
    function rotateEventTouch(ev){
       ev.preventDefault()
       currentX = ev.touches[0].clientX
  	   currentY = ev.touches[0].clientY
  	   var dX = currentX - previousX,
  	   dY = currentY - previousY;
  	   rotateData(2.0 * dX, 2.0 * dY)
  	   previousX = currentX;
       previousY = currentY;
    }
    
    function stopRotationTouch(ev){
      ev.preventDefault()
      canvas.removeEventListener('touchmove', rotateEventTouch)
  	 canvas.removeEventListener('touchend', stopRotationTouch)
  	 canvas.removeEventListener('touchcancel', stopRotationTouch)
    }
     
    function startRotatingTouch(ev){
       ev.preventDefault()
       previousX = ev.touches[0].clientX
  	   previousY = ev.touches[0].clientY
       canvas.addEventListener('touchmove', rotateEventTouch)
  	   canvas.addEventListener('touchend', stopRotationTouch)
  	   canvas.addEventListener('touchcancel', stopRotationTouch)
    }
  	
  	canvas.addEventListener('mousedown', startRotating)
     canvas.addEventListener('touchstart', startRotatingTouch)
     
     
gene_expression_div.innerHTML = ""
current_gene = geneselector.value;
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
    gene_raw_expression = []
    response.forEach(function(val){
      gene_raw_expression.push(parseFloat(val))
    })
    
    var vector = adaptiveMovingAverage(gene_raw_expression)
    vector = minMaxNormalization(vector)
    var gene_colors = valuesToColours(vector, 0, 1)
    var max_shown_value = 1
    plotGeneExpression(vector, max_shown_value, "Normalized gene expression")
  }
}
    queryString = "?gene_name=" + current_gene;
		request.open("GET", "fetch_gene_expression.php" + queryString, true)
		request.send(null)
  </script>
  <div style="float: clear;""><hr><span style="font-size:0.8em;">This data portal was created using the pseudotime_webportal tool (<a href="https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data">github link</a>) developed by Dorin-Mirel Popescu</span><hr></div>
</body>
</html>
"""

fetch_gene_expression = """
<?php 
	$gene_name  = $_GET["gene_name"];
	$gene_expression_file = './genes/' . $gene_name;
	$gene_expression_file = fopen($gene_expression_file, 'r');
	$gene_expression = fgets($gene_expression_file);
	echo $gene_expression;
?>
"""

import pandas as pd
import numpy as np

data = pd.read_csv(load_from, index_col = None)

# convert Colours to r, g, b values, then to floats < 1.0
def hexdec_to_1floats(hexdec):
    return np.array([int(hexdec[1:][i:(i+2)], 16) for i in (0, 2, 4)]) / 255.0

# get the limits fo pdt transitions
cell_labels  = data["Labels"]
Pdt          = data["Pseudotime"]
cell_colours = data["Colours"]
lab_centers  = []
cell_labels_unique = [l for l in np.unique(cell_labels)]
for l in cell_labels_unique:
    lab_centers.append(Pdt[cell_labels == l].median())
cell_labels_unique = [x for _, x in sorted(zip(lab_centers, cell_labels_unique))]
labels = cell_labels_unique
qqs = []
unique_colors = []
for l in cell_labels_unique:
    qqs.extend([q for q in np.percentile(Pdt[cell_labels == l].values, (5, 95))])
    unique_colors.append(cell_colours[cell_labels == l].values[0])
lab_transitions = []
for i in range(1, len(qqs) - 1, 2):
    lab_transitions.append((qqs[i] + qqs[i+1]) / 2)
lab_centers = []
for i in range(len(lab_transitions) + 1):
    if i == 0:
        first, last = 0, lab_transitions[i]
    elif i == len(lab_transitions):
        first, last = lab_transitions[i-1], 1
    else:
        first, last = lab_transitions[i-1], lab_transitions[i]
    lab_centers.append((first+last) / 2)
    
cell_labels_unique = ['"%s"' % l for l in cell_labels_unique]
cell_labels_unique = ",".join(cell_labels_unique);
cell_labels_unique = "var cell_labels_unique = [{d}]".format(d = cell_labels_unique)
lab_centers = sorted(lab_centers)
lab_centers = ",".join([str(l) for l in lab_centers])
lab_centers = "var lab_centers = [{d}]".format(d = lab_centers)
lab_transitions = ",".join([str(d) for d in lab_transitions])
lab_transitions = "lab_transitions = [{d}]".format(d = lab_transitions)
n_unique_colors = []
for n in unique_colors:
    rba = [str(i) for i in hexdec_to_1floats(n)]
    rba.append("0.2")
    rba = ",".join(rba)
    rba = "rgba({c})".format(c=rba)
    n_unique_colors.append(rba)
unique_colors = n_unique_colors
unique_colors = ['"{c}"'.format(c = c) for c in unique_colors]
unique_colors = ",".join(unique_colors)
unique_colors = "unique_colors = [{c}]".format(c = unique_colors)

# map Labels to colours
index_table = []
radio_commands = []
for index, label in enumerate(labels):
    indices = data.Labels == label
    indices = indices.values
    indices = np.where(indices)
    indices = ','.join([str(i) for i in indices[0]])
    indices = "[{indices}]".format(indices = indices)
    index_table.append("index_table['{label}'] = {indices}".format(label = label, indices = indices))
    colour = data.Colours[data.Labels == label].values[0]
    radio_command = "<div style='background-color:{colour}'><input type='checkbox' id='{label}' checked onchange='update3DPlot()' /><label for='{label}'>{label}: </label><br /></div>".format(colour = colour, label = label)
    radio_commands.append(radio_command)
index_table = ';\n  	'.join(index_table)
radio_commands = '\n					'.join(radio_commands)

# make data string
coordinates = data.values[:, 0:3].astype('float32')
# next few steps are compressing the data into a stadard cube centered at (0,0,0) and L = 200 
Xrange = np.percentile(coordinates[:, 0],  q = [1, 99]) * 1.2
Yrange = np.percentile(coordinates[:, 1],  q = [1, 99]) * 1.2
Zrange = np.percentile(coordinates[:, 2],  q = [1, 99]) * 1.2
center = np.tile(np.array([np.mean(Xrange), np.mean(Yrange), np.mean(Zrange)]), 
                 (coordinates.shape[0], 1))
coordinates = coordinates - center
Xrange = Xrange[1] - Xrange[0]
Yrange = Yrange[1] - Yrange[0]
Zrange = Zrange[1] - Zrange[0]
maxRange = max((Xrange, Yrange, Zrange))
ratio = 180.0 / maxRange
coordinates = coordinates * ratio

# next few steps the buffer data is created as string
colours     = data.values[:, 4]
buffer_data = []
for index in range(coordinates.shape[0]):
    coordinate = [str(i) for i in coordinates[index, :]]
    colour = [str(i) for i in hexdec_to_1floats(colours[index]).astype('float32')]
    vertex_data = coordinate + colour
    buffer_data.append(",".join(vertex_data))
buffer_data = ",".join(buffer_data)

pseudotime = data.values[:, 5]
pseudotime_buffer = []
for index in range(pseudotime.shape[0]):
    pseudotime_buffer.append(str(pseudotime[index, ]))
pseudotime_buffer = ",".join(pseudotime_buffer)

raw_expression = data.values[:, 6:]
gene_names     = data.columns[6:]
# make gene expression files here
gene_options = []
for index in range(gene_names.shape[0]):
    gene_name = gene_names[index]
    gene_expression = [str(val) for val in raw_expression[:, index]]
    gene_expression = ",".join(gene_expression)
    gene_name = gene_name.replace("/", "-")
    gene_expression_fpath = join(output_folder, "genes", gene_name)
    with open(gene_expression_fpath, "w") as gene_expression_fobj:
        gene_expression_fobj.write(gene_expression)
    gene_options.append("<option value='{gn}'>{gn}</option>".format(gn = gene_name));
gene_options = "".join(gene_options)
maxRawExpression = raw_expression.max()

template_str = template.replace('datahere', buffer_data)
template_str = template_str.replace('insert_lineage', insert_lineage)
template_str = template_str.replace('indiceshere', index_table)
template_str = template_str.replace('radiocommands', radio_commands)
template_str = template_str.replace('pseudotime_here', pseudotime_buffer)
template_str = template_str.replace('maxRawExpression_here', str(maxRawExpression))
template_str = template_str.replace('gene_options_here', gene_options)
template_str = template_str.replace('current_gene_here', "var current_gene = '{gn}'".format(gn = str(gene_names[0])))
template_str = template_str.replace('lab_centers_here', lab_centers)
template_str = template_str.replace('cell_labels_unique_here', cell_labels_unique)
template_str = template_str.replace('lab_transitions_here', lab_transitions)
template_str = template_str.replace('unique_colors_here', unique_colors)
with open(save_to, 'w') as result:
    result.write(template_str)
    
with open(fetch_php, "w") as fetch_php:
    fetch_php.write(fetch_gene_expression)
     