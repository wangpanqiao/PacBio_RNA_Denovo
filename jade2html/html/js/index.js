$(function(){
	if (navigator.userAgent.indexOf("MSIE 7.0")>0||navigator.userAgent.indexOf("MSIE 8.0")>0){
		$(".bs-sidebar").css("display","none");
		$("body").css("width","1024px");
	};
	$("html").find(".table").each(function(){
		$(this).find("tr:eq(0)").css("font-weight","bold").css("background-color","#333");
	});
	if($("h1").text()=="外显子分析报告"){
		var cols = $("table:eq(8) tr:eq(0) td").length;
		for(var i=18;i<26;i++){
			$("table:eq(8) tr:eq("+i+") td").attr("colspan",cols);
		}
	}
	for(var i=1;i<50;i++){
		$("#carousel"+i+"_t").flexslider({
			animation: "slide",
		    controlNav: false,
		    animationLoop: false,
		    slideshow: false,
		    itemWidth: 100,
		    itemMargin: 5,
		    asNavFor: "#carousel"+i
		});
		$("#carousel"+i).flexslider({
			animation: "slide",
		    controlNav: false,
		    animationLoop: false,
		    slideshow: false,
			sync: "#carousel"+i+"_t"
		});
	}
	$('a.lightbox').lightbox();
})
