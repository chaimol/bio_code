/**
 * Created by jian on 2016/9/30.
 */

$(document).ready(function() {
    /* 控制目录高度，当窗口高度改变时，目录高度也改变 */
    $(".catalog > nav").height( $(window).height()  - 100);
    $(window).resize(function(){
        $(".catalog > nav").height( $(window).height()  - 100);
    })

    /* 切换页面时记录当前高度，再次切换回来时滚动到此高度 */
    var tab_pane_offset={};
    $('#sidebar a[data-toggle="pill"]').each(function(){
        var tab_pane=$(this).attr('href');
        tab_pane_offset[tab_pane]=0;
    });
    $('#sidebar a[data-toggle="pill"]').on('hide.bs.tab', function (e) {
        var current=$(e.target).attr('href');
        var next=$(e.relatedTarget).attr('href');
        tab_pane_offset[current]=$(window).scrollTop();
        $('html, body').animate({scrollTop: tab_pane_offset[next] }, 0);
    })    

    //floatbar control
    $("#floatbar_goTop").click(function () {
        $('html, body').animate({scrollTop:0}, 200);
    })
    $(".float_btn").click(function () {
        var i=$('.float_btn').index($(this));
        if($('.float_tips').eq(i).hasClass('float_show')){
            $('.float_tips').eq(i).removeClass('float_show');
            $(".float_btn").eq(i).css('background-color','white');
            $(".float_btn span").eq(i).css('color','#0088cc');
            $(".float_btn p").eq(i).css('color','#0088cc');
            $(".float_btn p").eq(i).children("a").css('color','#0088cc');
        }else{
            $('.float_tips').removeClass('float_show');
            $(".float_btn").css('background-color','white');
            $(".float_btn span").css('color','#0088cc');
            $(".float_btn p").css('color','#0088cc');
            $(".float_btn p").children("a").css('color','#0088cc');
            $('.float_tips').eq(i).addClass('float_show');
            $(".float_btn").eq(i).css('background-color','#0088cc');
            $(".float_btn span").eq(i).css('color','white');
            $(".float_btn p").eq(i).css('color','white');
            $(".float_btn p").eq(i).children("a").css('color','white');
        }
    })
    $("#floatbar_remove").click(function () {
        $(".float_tips").css('display','none');
        $("#floatbar").css('display','none');
    })

  /*解析tables_js*/
  var tables_string=tables.toString().split('\n').slice(1,-1).join('\n') + '\n';
  var regex =/^\s*\[\s*([^\]]*)\s*\]\s*$/;
  var table_obj = {};
  var lines = tables_string.split(/\r\n|\r|\n/);
  var table_id;
  lines.forEach(function(line){
        //line=line.replace(/(^\s*)|(\s*$)/,'');
        if(regex.test(line)){
          var match = line.match(regex);
          table_obj[match[1]] = '';
          table_id = match[1];
        }else if(line.length != 0 ){
          table_obj[table_id] +=line+"\n";
        }
  });

  /*pdf首页*/
  var pdfFrontPage=table_obj["pdf"].split(/\n/);
  $("#project_name").html(pdfFrontPage[0]);
  $("#client").html(pdfFrontPage[1]);
  $("#date").html(pdfFrontPage[2]);
  delete table_obj["pdf"];

    /*有replace属性的元素按replace的值在tables.js里的json替换*/
    var replace_obj={};
    table_obj["replace_obj"].split(/\n/).map(function(line){ 
            if(line.length==0)return;
            var regex=/^\s*([\w\.\-\_]+)\s*=\s*(.*?)\s*$/;
            var match = line.match(regex);
            replace_obj[match[1]] = match[2].split(/,/);
    });
    delete table_obj["replace_obj"];
    $("[replace]").each(function(){
                var arr=replace_obj[$(this).attr("replace")];
                for (var i = 0; i < arr.length; i++) {
                    re = new RegExp($(this).attr("replace"),"gm");
                    $(this).before( $(this).prop("outerHTML").replace(re, arr[i]).replace(/replace=\"\w+\"/,''));
                }
                $(this).remove();
    })

    $('.imgbox').each(function(){
            $(this).find('a:first').tab('show');
    })

    /*为所有class 为table 的table元素创建表格*/
    $("[component_type=table] table").each(function(){
            if($(this).attr("id")){
                var rows= table_obj[$(this).attr("id")].split(/\n/);
                var tbody=$("<tbody/>");
                for (var i = 0; i < rows.length; i++) {
                    if(rows[i].length==0) continue;
                    var td = rows[i].split(/\t/);
                    var tr_obj=$("<tr/>");
                    for (var j = 0; j < td.length; j++) {
                        var td_obj = i==0 ? $("<th/>").html(td[j]) : $("<td/>").html(td[j]);
                        tr_obj.append(td_obj);
                    }
                        tbody.append(tr_obj);
                }
                $(this).append(tbody);
            }
    });

    /*Tips 点击跳转*/
    $('[data-toggle=tooltip]').click(function () {
        var $content=$($(this).attr('href')).parents('.content');
        var index=$(".content").index($content);
        $('#sidebar li:eq('+index+') a').tab('show');
        $('html, body').animate({scrollTop: $($(this).attr('href')).offset().top-50}, 400);
    });
});
//build orgchart_tree functions
function createNode(nodeData,level) {
    // construct the submodule of node
    var $nodeDiv = $('<div>').addClass('node');

    a_tips=nodeData["active"] != undefined ? '<a class="tips"  data-toggle="tooltip" title="点击跳转到对应章节。" href="#'+nodeData["href"]+'">' + nodeData["module"] + '</a>' : nodeData["module"];
    usable_class=nodeData["active"] != undefined ? 'usable':'unusable';
    if(nodeData["submodule"] != undefined){
        $nodeDiv.append('<div class="module '+usable_class+'">'+a_tips+'</div>')
        .append( '<div class="submodule">' + nodeData["submodule"] + '</div>' );
    }else{
        $nodeDiv.append('<div class="module_no_sub '+usable_class+'">'+a_tips+'</div>') ;
    }

    if (level > 0) {
      $($nodeDiv.children("div").get(0)).before('<i class="arrow_down"></i>');
    }
     
    return $nodeDiv;
  }
  // recursively build the tree
  function buildNode (nodeData, $appendTo, level) {
    var $table = $("<table cellpadding='0' cellspacing='0' border='0'/>");
    var $tbody = $("<tbody/>");

    // Construct the node
    var $nodeRow = $("<tr/>").addClass("node-cells");
    var $nodeCell = $("<td/>").addClass("node-cell").attr("colspan", 2);
    var $childNodes = nodeData.children;
    if ($childNodes && $childNodes.length > 1) {
      $nodeCell.attr("colspan", $childNodes.length * 2);
    }

    var $nodeDiv = createNode(nodeData, level);
    $nodeCell.append($nodeDiv);
    $nodeRow.append($nodeCell);
    $tbody.append($nodeRow);

    if ($childNodes && $childNodes.length > 0) {
      // recurse until leaves found (-1) or to the level specified
      var $childNodesRow;
      var $downLineRow = $("<tr/>");
      var $downLineCell = $("<td/>").attr("colspan", $childNodes.length * 2);
      $downLineRow.append($downLineCell);

      // draw the connecting line from the parent node to the horizontal line
      var $downLine = $("<div></div>").addClass("down");
      $downLineCell.append($downLine);
      $tbody.append($downLineRow);

      // draw the horizontal lines
      var $linesRow = $("<tr/>");
      $.each($childNodes, function() {
        var $left = $("<td>&nbsp;</td>").addClass("right top");
        var $right = $("<td>&nbsp;</td>").addClass("left top");
        $linesRow.append($left).append($right);
      });

      // horizontal line shouldn't extend beyond the first and last child branches
      $linesRow.find("td:first").removeClass("top").end().find("td:last").removeClass("top");
      $tbody.append($linesRow);

      $childNodesRow = $("<tr/>");
      $.each($childNodes, function() {
        var $td = $("<td class='node-container'/>");
        $td.attr("colspan", 2);
        // recurse through children lists and items

        buildNode(this, $td, level + 1);
        $childNodesRow.append($td);

      });
      $tbody.append($childNodesRow);
    }

    $table.append($tbody);
    $appendTo.append($table);

  };