var jade = require('jade'),
	fs = require("fs"),
	path = require("path");

var workdir = process.argv[2] || process.cwd();
var current_dir =process.cwd();
//workdir=path.relative(process.cwd(),workdir);
//console.log(workdir+'\n'+path.dirname(process.argv[1])+'\n');

var exec = require('child_process').exec,
    child;
var dir1 = path.dirname(process.argv[1])+"/public/";
var dir2 = path.dirname(process.argv[2])+"/";
var dir3 = dir2+"html";

//console.log(process.cwd());
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
    //console.log("!!!"+e);
})

// if(fs.existsSync(dir2))
// 	child = exec('rm -r '+dir2,function(error,stdout,strerr){
// 		console.log('rm -r '+dir2);
// 		if (error !== null) {
// 	      console.log('exec error: ' + error);
// 	    }
// 	});
child = exec('mkdir -p '+dir2+'/html'+' && cp -r '+dir1+'css '+dir3+' && cp -r '+dir1+'js '+dir3+' && cp -r '+dir1+'img '+dir3,
	function(error, stdout, stderr){
		// console.log('stdout: ' + stdout);
	 //    console.log('stderr: ' + stderr);
	 //    console.log('mkdir -p '+dir2+"/html");
	 //    console.log('cp -r '+dir1+"css "+dir2+'/html && cp -r '+dir1+'js '+dir3+' && cp -r'+dir1+'img '+dir3);
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

function Loadtxt(filename,R){
	var file_content = fs.readFileSync(filename.replace(/\\/g,"/"));
	var txt_obj = new Object();
	txt_obj.filename = filename;
	txt_obj.array = file_content.toString().split("\n");
	for(i in txt_obj.array) {
		if (txt_obj.array[i]=="") continue;
		txt_obj.array[i]=txt_obj.array[i].split(/[\t]+/);
	}
	R.push(txt_obj);
}

function getFiles(dir,result){
	var files = fs.readdirSync(dir);
	for(var i in files){
		if (!files.hasOwnProperty(i)) continue;
		var name = dir+"/"+files[i];
		if (fs.statSync(name).isDirectory()){
			getFiles(name,result);
		}else{
			if (name.search(/PacBio_rawdata.1.png/)!=-1){
				result.PacBio_rawdata.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/roi_stat.xls/)!=-1){
				Loadtxt(name,result.reads_of_insert_report);
			}
			else if (name.search(/cluster_stat.xls/)!=-1){
				Loadtxt(name,result.cluster_summary);
			}
			else if (name.search(/classify_stat.xls/)!=-1){
				Loadtxt(name,result.classify_summary_t);
			}
			else if (name.search(/Function_Annotation.stat.xls/)!=-1){
				Loadtxt(name,result.Function_Annotation);
			}
			else if (name.search(/KEGG_classification_count.txt/)!=-1){
				Loadtxt(name,result.kegg);
			}
			else if (name.search(/Total_mapped.stat.xls/)!=-1){
				Loadtxt(name,result.align_summary);
			}
			else if (name.search(/TranscriptExp.xls/)!=-1){
				Loadtxt(name,result.genes);
			}
			else if (name.search(/Kegg.gene2path.xls/)!=-1){
				Loadtxt(name,result.gene2path);
			}
			else if (name.search(/PacBio_sample.xls/)!=-1){
				Loadtxt(name,result.PacBio_sample);
			}
			else if (name.search(/group.xls/)!=-1){
				Loadtxt(name,result.group);
			}
			else if (name.search(/roi_readlength_hist.png/)!=-1){
				result.roi_readlength_hist.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/roi_accuracy_hist.png/)!=-1){
				result.roi_accuracy_hist.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/roi_npasses_hist.png/)!=-1){
				result.roi_npasses_histl4.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/classify_summary.png/)!=-1){
				result.classify_summary.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/ssr_density.png/)!=-1){
				result.ssr_density.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/Venn.png/)!=-1){
				result.venn.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/GO_classification_bar.png/)!=-1){
				result.go.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/UniIso.fa.KOG.class.png/)!=-1){
				result.kog.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/MA.png/)!=-1){
				result.MA.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/Volcano.png/)!=-1){
				result.Volcano.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/matrix.png/)!=-1){
				result.matrix.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/heatmap.png/)!=-1){
				result.heatmap.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/GO_Term.png/)!=-1){
				result.GO_Term.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/SigOfNodes.png/)!=-1){
				result.SigOfNodes.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/roi_length_density.png/)!=-1){
				result.roi_length_density.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/flnc_length_density.png/)!=-1){
				result.flnc_length_density.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/consensus_length_density.png/)!=-1){
				result.consensus_length_density.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/lncRNA_length_cout.png/)!=-1){
				result.lncRNA_length_cout.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/pathway_enrichment.png/)!=-1){
				result.pathway_enrichment.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/KEGG_classification.png/)!=-1){
				result.KEGG_classification.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/ko00010.png/)!=-1){
				result.KEGG_ko.push(name.replace(/\\/g,"/"));
			}
			else if (name.search(/_as.png/)!=-1){
				result.AS.push(name.replace(/\\/g,"/"));
			}
		//	else if (name.search(/.gif/)!=-1){
		//		result.PCA.push(name.replace(/\\/g,"/"));
		//	}
		//	else if (name.search(/down_and_down.png/)!=-1){
		//		result.profile.push(name.replace(/\\/g,"/"));
		//	}
		//	else if (name.search(/sampleClustering.png/)!=-1){
		//		result.sampleClustering.push(name.replace(/\\/g,"/"));
		//	}
		//	else if (name.search(/genes_clustering_dendrogram.png/)!=-1){
		//		result.genes_clustering_dendrogram.push(name.replace(/\\/g,"/"));
		//	}
		//	else if (name.search(/Module_trait_associations.png/)!=-1){
		//		result.Module_trait_associations.push(name.replace(/\\/g,"/"));
		//	}
		//	else if (name.search(/cytoscape_network\/.+.png/)!=-1){
		//		result.cytoscape_network.push(name.replace(/\\/g,"/"));
		//	}
			else if (name.search(/roi_length_density.png/)!=-1){
				result.roi_length_density.push(name.replace(/\\/g,"/"));
			}
		}
	}
}

var results=new Object();

results.roi_readlength_hist=[];
results.roi_accuracy_hist=[];
results.roi_npasses_histl4=[];
results.roi_length_density=[];

results.PacBio_rawdata=[];
results.cluster_summary=[];
results.ssr_density=[];
results.venn=[];
results.go=[];
results.kog=[];
results.MA=[];
results.Volcano=[];
results.heatmap=[];
results.matrix=[];
results.GO_Term=[];
results.SigOfNodes=[];
results.reads_of_insert_report=[];
results.classify_summary=[];
results.classify_summary_t=[];
results.Function_Annotation=[];
results.kegg=[];
results.align_summary=[];
results.genes=[];
results.gene2path=[];
results.PacBio_sample=[];
results.roi_length_density=[];
results.flnc_length_density=[];
results.consensus_length_density=[];
results.group=[];
results.lncRNA_length_cout=[];
results.KEGG_classification=[];
results.pathway_enrichment=[];
results.KEGG_ko=[];
results.PCA=[];
results.profile=[];
results.sampleClustering=[];
results.genes_clustering_dendrogram=[];
results.Module_trait_associations=[];
results.cytoscape_network=[];
results.AS=[];

chdir(path.dirname(workdir));
getFiles(path.basename(workdir),results);
fs.writeFile(dir2+"data.json", JSON.stringify(results), 'utf8', function(e){
	if(e)
		throw e;
	//fs.closeSync(fd);
});
//console.log(JSON.stringify(results.result,null, '\t'));
var config = require(current_dir+"/base_info.json");




var outfile=dir2;
var program_path=path.dirname(process.argv[1]);

jade.renderFile(path.join(program_path,'views/pac_rna_denovo.jade'), {filename:path.join(program_path,'views/pac_rna_denovo.jade'),
																  pretty:true,
																  pageTitle:'PacBio无参全长转录组测序分析报告',
																  data:JSON.stringify(results),
																  info:JSON.stringify(config)
															},
function (e, html) {
	if(e) throw e;
	fs.open(outfile+"report.html","w",0644,function(e,fd){
		console.log(outfile);
		if(e) throw e;
		fs.write(fd,html,0,'utf8',function(e){
			if(e) throw e;
			fs.closeSync(fd);
		})
	})
});

//jade.renderFile(path.join(program_path,'views/exon-info.jade'), {filename:path.join(program_path,'views/exon-info.jade'),
//																  pretty:true,
//																  pageTitle:'外显子流程结果说明文档'
//															},
//function (e, html) {
//	if(e) throw e;
//	fs.open(outfile+"ReadMe.html","w",0644,function(e,fd){
//		//console.log(outfile);
//		if(e) throw e;
//		fs.write(fd,html,0,'utf8',function(e){
//			if(e) throw e;
//			fs.closeSync(fd);
//		})
//	})
//});
chdir(current_dir);
