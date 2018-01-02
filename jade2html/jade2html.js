var jade = require('jade'),
	fs = require("fs"),
	path = require("path");

var workdir = process.argv[2] || process.cwd();
var current_dir =process.cwd();
workdir=path.relative(process.cwd(),workdir);
//console.log(workdir+'\n'+path.dirname(process.argv[1])+'\n');

var exec = require('child_process').exec,
    child,child1;
var dir1 = "/home/shenhc/nodeProject/jade2html/public/";
var dir2 = process.argv[3]||"/home/shenhc/report/";
var dir3 = dir2+"/html";

console.log(process.cwd());
//如果report存在则删除report
var rmdirSync = (function(){
    function iterator(url,dirs){
        var stat = fs.statSync(url);
        if(stat.isDirectory()){
            dirs.unshift(url);//收集目录
            inner(url,dirs);
        }else if(stat.isFile()){
            fs.unlinkSync(url);//直接删除文件
        }
    }
    function inner(path,dirs){
        var arr = fs.readdirSync(path);
        for(var i = 0, el ; el = arr[i++];){
            iterator(path+"/"+el,dirs);
        }
    }
    return function(dir,cb){
        cb = cb || function(){};
        var dirs = [];
        try{
            iterator(dir,dirs);
            for(var i = 0, el ; el = dirs[i++];){
                fs.rmdirSync(el);//一次性删除所有收集到的目录
            }
            cb()
        }catch(e){//如果文件或目录本来就不存在，fs.statSync会报错，不过我们还是当成没有异常发生
            e.code === "ENOENT" ? cb() : cb(e);
        }
    }
})();
rmdirSync(dir3,function(e){
    console.log("!!!"+e);    
})

// if(fs.existsSync(dir2))
// 	child = exec('rm -r '+dir2,function(error,stdout,strerr){
// 		console.log('rm -r '+dir2);
// 		if (error !== null) {
// 	      console.log('exec error: ' + error);
// 	    }
// 	});
child = exec('mkdir -p '+dir2+'/html'+' && cp -r '+dir1+'css '+dir3+' && cp -r '+dir1+'js '+dir3+' && cp -r '+dir1+'img '+dir3+' && cp -r '+workdir+' '+dir2,
	function(error, stdout, stderr){
		// console.log('stdout: ' + stdout);
	 //    console.log('stderr: ' + stderr);
	 console.log('mkdir -p '+dir2+"/html");
	 //    console.log('cp -r '+dir1+"css "+dir2+'/html && cp -r '+dir1+'js '+dir3+' && cp -r'+dir1+'img '+dir3);
	 console.log(workdir+' '+dir2);
	    if (error !== null) {
	      console.log('exec error: ' + error);
	    }
	}
);
function chdir(dir){
	try {
		process.chdir(dir);
		//console.log('New directory: ' + process.cwd());
	}
	catch (err) {
		console.log('chdir: ' + err);
	}
}

//process.argv[2]的上一级目录
// function whoisf(w){
// 	var last = w.lastIndexOf('\\');
// 	console.log(last);
// 	return w.substr(0,last);	
// }
// var fwork = whoisf(workdir);
// var fwork = path.normalize(workdir+"\\..");
// //craete /html/css&&/js
// function mkdirSync(url,mode,cb){
//     var arr = url.split("/");
//     mode = mode || 0755;
//     cb = cb || function(){};
//     if(arr[0]==="."){//处理 ./aaa
//         arr.shift();
//     }
//     if(arr[0] == ".."){//处理 ../ddd/d
//         arr.splice(0,2,arr[0]+"/"+arr[1])
//     }
//     function inner(cur){
//         if(!fs.existsSync(cur)){//不存在就创建一个
//             fs.mkdirSync(cur, mode)
//         }
//         if(arr.length){
//             inner(cur + "/"+arr.shift());
//         }else{
//             cb();
//         }
//     }
//     arr.length && inner(arr.shift());
// }
// mkdirSync(fwork+"\\html", 755,function(e){
// 	if(e){
//         console.log("创建"+fwork+'\\html 失败。');
//     }
// });
// mkdirSync(fwork+"\\html\\css", 755,function(e){
// 	if(e)
// 		console.log("创建"+fwork+'\\html\\css 失败。');
// });
// mkdirSync(fwork+"\\html\\js", 755,function(e){
// 	if(e)
// 		console.log("创建"+fwork+'\\html\\js 失败。');
// });
// mkdirSync(fwork+"\\html\\img", 755,function(e){
// 	if(e)
// 		console.log("创建"+fwork+'\\html\\img 失败。');
// });
// //copy files
// function copy(filename){
// 	var readStream = fs.createReadStream(process.cwd()+"\\public\\"+filename);
// 	var writeStream = fs.createWriteStream(fwork+"\\html\\"+filename);
// 	console.log(process.cwd()+"\\public\\"+filename);
// 	console.log(fwork+"\\html\\"+filename);
// 	readStream.pipe(writeStream);
// 	readStream.on('end', function () {
// 	 console.log('copy end');
// 	});
// 	readStream.on('error', function () {
// 	 console.log('copy error');
// 	});
// }
// var files = ["js\\jquery-1.10.2.js","js\\bootstrap.js","js\\index.js","css\\bootstrap.css","css\\index.css","img\\RNA-seq_bioinfo_pipeline.png","img\\RNA-seq_library_construction.png"];
// for (i in files){
// 	copy(files[i]);
// };
function Loadtxt(filename,R){
	var file_content = fs.readFileSync(filename.replace(/\\/g,"/"));
	var array = file_content.toString().split("\n");
	for(i in array) {
		if (array[i]=="") continue;
		array[i]=array[i].split(/[\t|,]+/);
		R.push(array[i]);
	}
}

function getFiles(dir,result){
	var files = fs.readdirSync(dir);
	for(var i in files){
		if (!files.hasOwnProperty(i)) continue;
		var name = dir+'/'+files[i];		
		if (fs.statSync(name).isDirectory()){
			getFiles(name,result);
		}else{
			if (name.search(/alignment_stat.csv/)!=-1){
				Loadtxt(name,result.alignment_stat);
			}
			else if (name.search(/fastqc_report.html/)!=-1){
				result.fastqc_report.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/saturation_covered_genes.png/)!=-1){
				result.saturation_covered_genes.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/saturation_RPKM.png/)!=-1){
				result.saturation_RPKM.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/RPKM_range.png/)!=-1){
				result.RPKM_range.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/randCheck_mRNA.png/)!=-1){
				result.randCheck_mRNA.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/randCheck_gene.png/)!=-1){
				result.randCheck_gene.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/region_Pie_Chart.png/i)!=-1){
				result.region_Pie_Chart.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/sense_antisense_Pie_Chart.png/)!=-1){
				result.sense_antisense_Pie_Chart.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/quals.png/)!=-1){
				result.quals.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/quals2.png/)!=-1){
				result.quals2.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/quals3.png/)!=-1){
				result.quals3.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/quals-hm.png/)!=-1){
				result.quals_hm.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/acgt-cycles.png/)!=-1){
				result.acgt_cycles.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/gc-content.png/)!=-1){
				result.gc_content.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/gc-depth.png/)!=-1){
				result.gc_depth.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/insert-size.png/)!=-1){
				result.insert_size.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/coverage.png/)!=-1){
				result.coverage.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/mism-per-cycle.png/)!=-1){
				result.mism_per_cycle.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/indel-dist.png/)!=-1){
				result.indel_dist.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/indel-cycles.png/)!=-1){
				result.indel_cycles.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_DEG.png/)!=-1){
				result.DEG.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_[0-9]+_Up.txt/)!=-1){
				Loadtxt(name,result.DEG_UpList);
			}
			else if (name.search(/_[0-9]+_Down.txt/)!=-1){
				Loadtxt(name,result.DEG_DownList);
			}
			else if (name.search(/_reads_linear.png/)!=-1){
				result.reads_linear.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_GO_Term.png/)!=-1){
				result.GO_Term.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_SigOfNodes.png/)!=-1){
				result.SigOfNodes.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_GenTable.txt/)!=-1){
				Loadtxt(name,result.GenTable);
			}
			else if (name.search(/venn.png/)!=-1){
				result.venn.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_heatmap.png/)!=-1){
				result.heatmap.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/pathway_enrichment.png/)!=-1){
				result.pathway.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_ctrl_Down_pathway_enrichment.xls/)!=-1){
				Loadtxt(name,result.Down_pathway);
			}
			else if (name.search(/_ctrl_Up_pathway_enrichment.xls/)!=-1){
				Loadtxt(name,result.Up_pathway);
			}
			else if (name.search(/W_Gene_Regions_Info.png/)!=-1){
				result.Regions_Info.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/New_Gene_predict.gtf/)!=-1){
				Loadtxt(name,result.raw);
			}
		}
	}
}

var results=new Object();
results.alignment_stat=[];
results.fastqc_report=[];
results.saturation_covered_genes=[];
results.saturation_RPKM=[];
results.RPKM_range=[];
results.randCheck_mRNA=[];
results.randCheck_gene=[];
results.region_Pie_Chart=[];
results.sense_antisense_Pie_Chart=[];
results.quals=[];
results.quals2=[];
results.quals3=[];
results.quals_hm=[];
results.acgt_cycles=[];
results.gc_content=[];
results.gc_depth=[];
results.insert_size=[];
results.coverage=[];
results.venn=[];
results.mism_per_cycle=[];
results.indel_dist=[];
results.indel_cycles=[];
results.DEG=[];
results.DEG_UpList=[];
results.DEG_DownList=[];
results.reads_linear=[];
results.GO_Term=[];
results.GenTable=[]
results.SigOfNodes=[];
results.heatmap=[];
results.pathway=[];
results.Down_pathway=[]
results.Up_pathway=[];
results.Regions_Info=[];
results.raw=[];
chdir(path.dirname(workdir));
getFiles(path.basename(workdir),results);
// fs.writeFile(dir2+"data.json", JSON.stringify(results), 'utf8', function(e){
// 	if(e)
// 		throw e;
// 	//fs.closeSync(fd);
// });
//console.log(JSON.stringify(results));


var outfile=process.argv[3] || "/home/shenhc/report/";
var program_path=path.dirname(process.argv[1]);

//var config = require(path.relative(current_dir,workdir+"/base_info.json"));
var config = require(process.argv[2]+"/base_info.json");
//var config = require("/home/shenhc/nodeProject/jade2html/base_info.json");
jade.renderFile(path.join(program_path,'views/index.jade'), {filename:path.join(program_path,'views/index.jade'),
																  pretty:true,
																  pageTitle:'RNA-Seq数据分析报告',
																  data:JSON.stringify(results),
																  info:JSON.stringify(config)
															},
function (e, html) {
	if(e) throw e;
	fs.open(outfile+"/index.html","w",0644,function(e,fd){
		//console.log(outfile);
		//console.log(config.species);
		if(e) throw e;
		fs.write(fd,html,0,'utf8',function(e){
			if(e) throw e;
			fs.closeSync(fd);
		})
	})
})
var outfile_info=path.dirname(process.argv[1])+"/";

jade.renderFile(path.join(program_path,'views/rna-seq-info.jade'), {filename:path.join(program_path,'views/rna-seq-info.jade'),
																  pretty:true,
																  pageTitle:''
															},
function (e, html) {
	if(e) throw e;
	fs.open(outfile+"/ReadMe.html","w",0644,function(e,fd){
		//console.log(outfile_info);
		if(e) throw e;
		fs.write(fd,html,0,'utf8',function(e){
			if(e) throw e;
			fs.closeSync(fd);
		})
	})
});
chdir(current_dir);
//craete /html/css&&/js
// function mkdirSync(url,mode,cb){
//     var arr = url.split("/");
//     mode = mode || 0755;
//     cb = cb || function(){};
//     if(arr[0]==="."){//处理 ./aaa
//         arr.shift();
//     }
//     if(arr[0] == ".."){//处理 ../ddd/d
//         arr.splice(0,2,arr[0]+"/"+arr[1])
//     }
//     function inner(cur){
//         if(!fs.existsSync(cur)){//不存在就创建一个
//             fs.mkdirSync(cur, mode)
//         }
//         if(arr.length){
//             inner(cur + "/"+arr.shift());
//         }else{
//             cb();
//         }
//     }
//     arr.length && inner(arr.shift());
// }
// var fwork = path.normalize(workdir+"/..");
// mkdirSync(fwork+"/html", 755,function(e){
// 	if(e){
//         console.log("创建"+fwork+'/html 失败。');
//     }
// });

// child = exec('mkdir ../report && cp -r '+path.join(program_path,'public')+' '+outfile+' '+workdir+' report',
//   function (error, stdout, stderr) {
//   	console.log(path.join(program_path,'public')+' '+outfile+' '+workdir+' report');
//    console.log('stdout: ' + stdout);
//    console.log('stderr: ' + stderr);
//     if (error !== null) {
//       console.log('exec error: ' + error);
//     }
// });