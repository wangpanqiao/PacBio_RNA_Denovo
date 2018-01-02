var jade = require('jade'),
	fs = require("fs"),
	path = require("path");

var workdir = process.argv[2] || process.cwd();
workdir=path.relative(process.cwd(),workdir);
console.log(workdir+'\n'+path.dirname(process.argv[1])+'\n');

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
			else if (name.search(/region_Pie_Chart.png/)!=-1){
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
getFiles(workdir,results);
//console.log(JSON.stringify(results));

var outfile=process.argv[3] || "test.htm";
var program_path=path.dirname(process.argv[1]);
jade.renderFile(path.join(program_path,'views/index_v2.jade'), {filename:path.join(program_path,'views/index_v2.jade'),
																  pretty:true,
																  pageTitle:'RNA-Seq数据分析报告',
																  data:JSON.stringify(results)
															},
function (e, html) {
	if(e) throw e;
	fs.open(outfile,"w",0644,function(e,fd){
		if(e) throw e;
		fs.write(fd,html,0,'utf8',function(e){
			if(e) throw e;
			fs.closeSync(fd);
		})
	})
})

var exec = require('child_process').exec,
    child;

child = exec('mkdir report && cp -r '+path.join(program_path,'public')+' '+outfile+' '+workdir+' report'+' && rm '+outfile,
  function (error, stdout, stderr) {
//    console.log('stdout: ' + stdout);
//    console.log('stderr: ' + stderr);
    if (error !== null) {
      console.log('exec error: ' + error);
    }
});