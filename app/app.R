library(shiny)
library(shinythemes)
library(dplyr)
library(DT)
library(car)
library(ggplot2)
library(Matrix)
library(foreach)
library(doParallel)
library(grid)
library(webshot)
library(cowplot)
library(Cairo)
library(shinyjs)
library(openxlsx) # 엑셀 파일 생성
library(zip) # ZIP 파일 생성
library(igraph) # 그래프 생성 및 Leiden 클러스터링
library(visNetwork) # 인터랙티브 네트워크 그래프
library(RANN) # KD-tree nearest neighbor
library(RcppHNSW) # HNSW
library(matrixStats)
library(viridis)

## [1] ▶
example_data <- local({
  files <- c(expr   = "test_data/count_mat_p16.rda",
             coords = "test_data/coordinates_p16.rda",
             labels = "test_data/clusters_p16.rds")
  if (!all(file.exists(files))) {
    warning("Example data not found: ", paste(files[!file.exists(files)], collapse = ", "))
    return(NULL)
  }
  
  env <- new.env()
  load(files[["expr"]],   envir = env)   # -> env$count_mat_p16
  load(files[["coords"]], envir = env)   # -> env$coordinates_p16
  
  ## dgCMatrix -> dense matrix (앱의 CSV 업로드 경로와 동일한 형태: as.matrix(expr))
  expr <- as.matrix(env$count_mat_p16)
  storage.mode(expr) <- "double"
  
  ## 좌표: 앱이 기대하는 Spot / X / Y 형식으로 변환, expr 열 순서에 맞춰 정렬
  co <- env$coordinates_p16[colnames(expr), , drop = FALSE]
  coords <- data.frame(Spot = colnames(expr),
                       X    = co$imagerow,
                       Y    = co$imagecol,
                       stringsAsFactors = FALSE)
  
  ## BANYAN 라벨: named vector (names = barcode) 를 expr 열 순서로 재정렬
  ## unname() 필수 — 이름이 남아 있으면 sendCustomMessage 가 JSON 배열이 아닌 객체로 보냄
  labels <- unname(readRDS(files[["labels"]])[colnames(expr)])
  
  stopifnot(nrow(coords) == ncol(expr),
            !anyNA(coords$X), !anyNA(coords$Y), !anyNA(labels))
  
  list(expr = expr, coords = coords, labels = labels,
       description = sprintf("%s genes x %s spots, %d tissue domains",
                             format(nrow(expr), big.mark = ","),
                             format(ncol(expr), big.mark = ","),
                             length(unique(labels))))
})


# 커스텀 CSS
custom_css <- "
/* 기본 폰트 및 스타일 */
body, h3, h4, p {
  font-family: 'Arial', sans-serif;
}

/* 버튼 스타일 */
.btn-primary, .btn-danger, .btn-warning {
  transition: all 0.3s ease;
  border-radius: 5px;
}

.btn-primary:hover, .btn-danger:hover, .btn-warning:hover {
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
}

/* 툴팁 스타일 */
.custom-tooltip {
  position: absolute;
  background-color: #333;
  color: white;
  padding: 5px 10px;
  border-radius: 3px;
  font-size: 12px;
  z-index: 2000; /* 캔버스 위에 표시되도록 z-index 증가 */
  pointer-events: none;
  white-space: nowrap;
}

/* 사이드바 패널 */
.sidebar-panel {
  background-color: #f8f9fa;
  padding: 15px;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  overflow-y: auto;
}

/* 탭 콘텐츠 */
.tab-content {
  padding: 20px;
}

.tab-panel-title {
  font-size: 16px;
  color: #2980b9;
  font-weight: bold;
}

/* 캔버스 컨테이너 */
.canvas-container {
  border: 1px solid #ddd;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  display: block;
}

/* 조정 행 */
.adjustment-row {
  margin-bottom: 45px;
}

/* 숫자 입력 컨테이너 */
.numeric-input-container {
  margin-bottom: 10px;
}

/* 현재 클러스터 텍스트 */
.current-cluster-text {
  font-size: 14px;
  font-weight: 700;
  color: #2c3e50;
  margin-top: 10px;
}

/* Test Results 컨테이너 */
.test-results-container {
  margin-top: 20px;
}

/* 테이블 컨테이너 */
.table-container {
  background-color: #f9f9f9;
  padding: 15px;
  border-radius: 5px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
}

/* DataTable 스타일 */
.dataTable td {
  padding: 8px;
  font-size: 13px;
}

.dataTable th {
  background-color: #2980b9;
  color: white;
  font-weight: bold;
  padding: 8px;
}

/* 데이터 없음 메시지 */
.no-data {
  text-align: center;
  color: #7f8c8d;
  padding: 20px;
}

.dataTable th {
  background-color: #2980b9;
  color: white;
  font-weight: bold;
  padding: 6px;
}

/* 페이지 네비게이션 스타일 */
.dataTables_paginate .paginate_button {
  padding: 5px 10px;
  margin: 0 2px;
  border-radius: 3px;
  background-color: #f8f9fa;
  color: #2980b9;
}

.dataTables_paginate .paginate_button:hover {
  background-color: #2980b9;
  color: white;
}

.dataTables_paginate .paginate_button.current {
  background-color: #2980b9;
  color: white;
  font-weight: bold;
}

/* 다운로드 버튼과 테이블 간격 조정 */
.download-btn-container {
  margin-top: 15px;
  margin-bottom: 10px;
}

/* 테이블과 페이지네이션 간격 조정 */
.dataTables_wrapper .bottom {
  margin-top: 10px;
}

/* 카드 스타일 */
.gene-card {
  background-color: #ffffff;
  border: 1px solid #e0e0e0;
  border-radius: 6px;
  padding: 10px;
  margin-bottom: 8px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
  transition: transform 0.2s ease, box-shadow 0.2s ease;
  cursor: pointer;
}

.gene-card:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
}

.gene-card-header {
  font-size: 13px;
  font-weight: bold;
  color: #2c3e50;
  margin-bottom: 5px;
}

.gene-card-content {
  font-size: 11px;
  color: #7f8c8d;
}

.gene-card-content span {
  display: inline-block;
  margin-right: 10px;
}

/* 검색 바 스타일 */
.search-bar {
  width: 100%;
  padding: 8px;
  margin-bottom: 10px;
  border: 1px solid #ddd;
  border-radius: 4px;
  font-size: 12px;
  box-sizing: border-box;
}

/* 카드 리스트 컨테이너 */
.cards-container {
  max-height: 440px;
  overflow-y: auto;
  padding-right: 5px;
}

/* 스크롤바 스타일링 */
.cards-container::-webkit-scrollbar {
  width: 6px;
}

.cards-container::-webkit-scrollbar-track {
  background: #f1f1f1;
}

.cards-container::-webkit-scrollbar-thumb {
  background: #2980b9;
  border-radius: 3px;
}

.cards-container::-webkit-scrollbar-thumb:hover {
  background: #1f6391;
}

/* 클릭 가능한 유전자 이름 스타일 */
.clickable-gene {
  color: #1e90ff; /* 밝은 파란색 */
  font-weight: bold; /* 굵은 글씨로 강조 */
  text-decoration: none; /* 밑줄 제거 */
  background-color: #e6f3ff; /* 연한 배경색 */
  padding: 4px 8px; /* 패딩으로 버튼 느낌 */
  border-radius: 4px; /* 둥근 테두리 */
  border: 1px solid #1e90ff; /* 테두리 추가 */
  transition: all 0.2s ease; /* 부드러운 전환 효과 */
  display: inline-flex; /* 아이콘과 텍스트 정렬 */
  align-items: center; /* 세로 중앙 정렬 */
  cursor: pointer;
  position: relative; /* 툴팁 위치를 위해 */
}

/* 호버 효과 강화 */
.clickable-gene:hover {
  color: #ffffff; /* 흰색 텍스트 */
  background-color: #2980b9; /* 진한 파란색 배경 */
  border-color: #2980b9; /* 테두리 색상 일치 */
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2); /* 그림자 추가 */
}

/* 외부 링크 아이콘 */
.clickable-gene::after {
  content: '\\f35d'; /* Font Awesome 링크 아이콘 */
  font-family: 'Font Awesome 6 Free';
  font-weight: 900;
  font-size: 12px;
  margin-left: 5px;
  color: inherit;
}

/* 툴팁 스타일 */
.clickable-gene-tooltip {
  visibility: hidden;
  background-color: #333;
  color: white;
  text-align: center;
  padding: 5px 10px;
  border-radius: 3px;
  font-size: 12px;
  position: absolute;
  z-index: 1000;
  bottom: 125%; /* Gene 이름 위에 표시 */
  left: 50%;
  transform: translateX(-50%);
  white-space: nowrap;
  opacity: 0;
  transition: opacity 0.2s;
}

.clickable-gene:hover .clickable-gene-tooltip {
  visibility: visible;
  opacity: 1;
}

/* 유전자 분석 헤더 스타일 */
.gene-analysis-container h3 .clickable-gene {
  font-size: 18px;
  margin-bottom: 10px;
}

/* Tissue Analysis 컨테이너 */
.tissue-analysis-container {
  margin-top: 20px;
}

/* Friedman 결과 박스 */
.friedman-result-box {
  background-color: #f8f9fa;
  padding: 15px;
  border-radius: 5px;
  margin-top: 10px;
  border-left: 4px solid #2980b9;
}

/* 클러스터 카드 특별 스타일 */
.cluster-card-selected {
  background-color: #e8f4f8 !important;
  font-weight: bold !important;
  box-shadow: 0 2px 8px rgba(41, 128, 185, 0.3) !important;
  transform: translateY(-2px);
}

.cluster-card-current {
  background-color: #fff3cd !important;
  border-color: #ffeaa7 !important;
}

/* 카드 호버 효과 강화 */
#cluster-status .gene-card:hover {
  transform: translateY(-3px);
  box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
  border-left-width: 6px;
}

"

## [1] ▶▶
source("help_content.R", local = TRUE, encoding = "UTF-8")

# UI 정의
ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  tags$head(
    tags$style(HTML(custom_css)),
    tags$style(HTML(help_css)), ## [2] ▶▶
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"),
    tags$script(HTML('
      $(document).ready(function() {
        // 툴팁 요소 생성
        if ($("#customTooltip").length === 0) {
          $("body").append("<div id=\'customTooltip\' class=\'custom-tooltip\' style=\'display: none;\'></div>");
        }
        var $tooltip = $("#customTooltip");
        
        // 메시지 저장용 변수
        var tooltipMessages = {
          "test_button": "",
          "gene_button": "",
          "new_cluster_btn": "",
          "run_spatial_analysis": ""
        };
        
        // Shiny로부터 메시지 수신
        Shiny.addCustomMessageHandler("updateTooltipMessage", function(message) {
          tooltipMessages[message.buttonId] = message.text;
        });

        // 비활성화된 버튼에 마우스 이벤트 추가
        $(document).on("mouseenter", ".btn[disabled]", function(e) {
          var buttonId = $(this).attr("id");
          var message = tooltipMessages[buttonId] || "Unknown error";
          $tooltip.text(message);
          $tooltip.show();
          updateTooltipPosition(e);
        });

        $(document).on("mousemove", ".btn[disabled]", function(e) {
          updateTooltipPosition(e);
        });

        $(document).on("mouseleave", ".btn[disabled]", function() {
          $tooltip.hide();
        });

        function updateTooltipPosition(e) {
          var x = e.pageX + 10; // 커서 오른쪽으로 10px
          var y = e.pageY + 10; // 커서 아래로 10px
          $tooltip.css({ top: y + "px", left: x + "px" });
        }
      });
    '))
  ),
  
  titlePanel(
    div(
      strong("VeteranST", 
             style = "font-size: 24px; color: #ffffff; font-style: italic; font-weight: 500; font-family: 'Impact', 'Arial Black', sans-serif; letter-spacing: 4px;"),
      style = "text-align: center; padding: 15px; background-color: #2980b9; border-radius: 8px; margin-bottom: 20px;"
    )
  ),
  tabsetPanel(
    id = "main_tabs", ## [3] ▶▶
    type = "pills",
    tabPanel(
      tags$span(icon("upload"), " Data Upload", class = "tab-panel-title"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          step_header("STEP 1. Load Data", "step1"), ## [4] ▶▶
          fileInput("expr_file", "Expression Matrix(.csv)", accept = ".csv", width = "100%"),
          fileInput("coords_file", "Spot Coordinates(.csv)", accept = ".csv", width = "100%"),
          uiOutput("cluster_file_ui"),
          ## [2] ▶
          tags$hr(style = "margin: 18px 0 12px 0;"),
          tags$h4("Or use the example data",
                  style = "color: #2980b9; font-size: 15px; margin-bottom: 8px;"),
          actionButton("load_example_btn",
                       tags$span(icon("flask"), " Load Example Data"),
                       class = "btn-success btn-sm", width = "100%"),
          tags$p("Melanoma brain metastasis (P16; Sudmeier et al., 2022), 10x Visium.",
                 tags$br(),
                 "16,074 genes \u00d7 1,650 spots with BANYAN tissue domains.",
                 style = "font-size: 11px; color: #7f8c8d; margin-top: 8px; margin-bottom: 0;")
        ),
        mainPanel(
          width = 9,
          div(
            style = "display: flex; flex-direction: column; align-items: flex-start;",
            div(
              style = "position: relative; width: 500px; margin-left: 10px;",
              h3("Tissue Spots", style = "color: #2c3e50; font-size: 20px; position: absolute; top: -30px; left: 50%; transform: translateX(-50%); margin: 0;"),
              tags$canvas(id = "dataCanvas", width = "500", height = "500", class = "canvas-container")
            )
          )
        )
      )
    ),
    tabPanel(
      tags$span(icon("paint-brush"), " Tissue Domain Assignment", class = "tab-panel-title"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          step_header("STEP 2. Assign Tissue Domain", "step2"), ## [4] ▶▶
          uiOutput("new_cluster_btn_ui", style = "margin-bottom: 10px;"),
          div(style = "margin-top: 15px;"),
          uiOutput("erase_btn_ui", style = "margin-bottom: 10px;"),
          sliderInput("brush_size", "Brush/Eraser Size", min = 1, max = 50, value = 20, step = 5, width = "100%", ticks = FALSE),
          div(style = "margin-top: 15px;"),
          actionButton("auto_cluster_btn", tags$span(icon("magic"), " Set the rest as a new domain"), class = "btn-primary btn-sm", width = "100%", ticks = FALSE),
          div(style = "margin-top: 15px;"),
          actionButton("reset_cluster_btn", tags$span(icon("undo"), " Initialize"), class = "btn-danger btn-sm", width = "100%"),
          div(style = "margin-top: 30px;"),
          downloadButton("download_clusters", "Download Domain Assignments", class = "btn-success btn-sm", style = "width: 100%;")
        ),
        mainPanel(
          width = 9,
          div(
            style = "display: flex; flex-direction: column; align-items: flex-start;",
            div(
              style = "position: relative; width: 800px; margin-left: 10px;",
              h3("Tissue Domain Assignment", style = "color: #2c3e50; font-size: 20px; position: absolute; top: -30px; left: 250px; transform: translateX(-50%); margin: 0;"),
              div(
                style = "display: flex; align-items: flex-start;",
                tags$canvas(id = "spotCanvas", width = "500", height = "500", class = "canvas-container"),
                div(
                  style = "margin-left: 20px; width: 300px;",
                  uiOutput("cluster_status_section")
                )
              )
            )
          )
        )
      ),
      tags$head(
        tags$script(HTML('
          // 툴팁 초기화
          $(document).ready(function() {
            if ($("#customTooltip").length === 0) {
              $("body").append("<div id=\'customTooltip\' class=\'custom-tooltip\' style=\'display: none;\'></div>");
            }
          });

          var spotData = [];
          var brushSize = 20;
          var currentCluster = 0;
          var clusterInfo = {};
          var tempSelectedSpots = new Set();
          var eraseMode = false;
          var highlightTimeout = null;
          
            function getClusterColor(clusterIndex) {
              if (clusterIndex === 0) return "black";
              const hue = (clusterIndex * 137.5) % 360;
              return `hsl(${hue}, 70%, 50%)`;
            }

          
          function renderCanvas(canvasId, clear = true, showClusters = true, clusterData = null, highlightedCluster = null) {
            var canvas = document.getElementById(canvasId);
            var ctx = canvas.getContext("2d");
            if (clear) {
              ctx.clearRect(0, 0, canvas.width, canvas.height);
              drawAllSpots(canvasId, showClusters, clusterData, highlightedCluster);
            }
          }
          
          function drawAllSpots(canvasId, showClusters, clusterData, highlightedCluster) {
            var canvas = document.getElementById(canvasId);
            var ctx = canvas.getContext("2d");
            spotData.forEach(function(spot, index) {
              var cluster = (canvasId === "dataCanvas" && clusterData && clusterData[index] !== undefined) 
                ? clusterData[index] 
                : (spot.cluster || 0);
    
              ctx.beginPath();
              ctx.arc(spot.X, spot.Y, 3, 0, 2 * Math.PI);
    
              if (highlightedCluster !== null) {
                // 문자열과 숫자 비교를 위해 둘 다 숫자로 변환하여 비교
                var clusterNum = Number(cluster);
                var highlightedNum = Number(highlightedCluster);
      
                if (clusterNum === highlightedNum && clusterNum !== 0) {
                  ctx.fillStyle = getClusterColor(clusterNum);
                  ctx.fill();
                } else if (clusterNum !== 0) {
                  ctx.fillStyle = getClusterColor(clusterNum).replace(")", ", 0.4)").replace("hsl", "hsla");
                  ctx.fill();
                } else {
                  ctx.strokeStyle = "black";
                  ctx.lineWidth = 1;
                  ctx.stroke();
                }
              } else if (showClusters && cluster !== 0) {
                ctx.fillStyle = getClusterColor(cluster);
                ctx.fill();
              } else {
                ctx.strokeStyle = "black";
                ctx.lineWidth = 1;
                ctx.stroke();
              }
            });
          }
          
          Shiny.addCustomMessageHandler("highlightCluster", function(data) {
  
            // 기존 타이머가 있으면 취소
            if (highlightTimeout) {
              clearTimeout(highlightTimeout);
              highlightTimeout = null;
            }
  
            spotData = JSON.parse(JSON.stringify(data.spots));
            currentCluster = data.currentCluster;
            if (data.brushSize) brushSize = data.brushSize;
            clusterInfo = data.clusterInfo || {};
            
  
            // 하이라이트 렌더링
            renderCanvas("spotCanvas", true, true, null, data.highlightedCluster);
  
            // 1.5초 후 하이라이트 해제
            if (data.highlightedCluster !== null) {
              highlightTimeout = setTimeout(function() {
                console.log("Highlight timeout - removing highlight");
                renderCanvas("spotCanvas", true, true, null, null);
                highlightTimeout = null;
              }, 1500);
            }
          });
          
          
          var previousCluster = 0;

          Shiny.addCustomMessageHandler("initDataCanvas", function(data) {
            spotData = JSON.parse(JSON.stringify(data.spots));
            renderCanvas("dataCanvas", true, true, data.uploadedClusters);
          });
          
          Shiny.addCustomMessageHandler("initSpotCanvas", function(data) {
            spotData = JSON.parse(JSON.stringify(data.spots));
            brushSize = data.brushSize;
            currentCluster = data.currentCluster;
            clusterInfo = data.clusterInfo || {};
            tempSelectedSpots = new Set();
            previousCluster = currentCluster;
            eraseMode = false;
            renderCanvas("spotCanvas", true, true);

            var canvas = document.getElementById("spotCanvas");
            var ctx = canvas.getContext("2d");
            var drawing = false;
            var lastDrawnSpots = new Set();
            var $tooltip = $("#customTooltip");

            function selectSpots(cx, cy) {
              var newSpots = new Set();
              spotData.forEach(function(spot, index) {
                var x = spot.X;
                var y = spot.Y;
                var dist = Math.sqrt((x - cx) ** 2 + (y - cy) ** 2);
                if (dist <= brushSize) {
                  if (!eraseMode && spot.cluster === 0) {
                    newSpots.add(index);
                  } else if (eraseMode && spot.cluster === currentCluster) {
                    newSpots.add(index);
                  }
                }
              });
              newSpots.forEach(function(index) {
                tempSelectedSpots.add(index);
              });
              return newSpots;
            }

            canvas.onmousedown = function(e) {
              if (currentCluster === 0) return;
              drawing = true;
              var rect = canvas.getBoundingClientRect();
              var x = e.clientX - rect.left;
              var y = e.clientY - rect.top;
              selectSpots(x, y);
            };

            canvas.onmousemove = function(e) {
              var rect = canvas.getBoundingClientRect();
              var x = e.clientX - rect.left;
              var y = e.clientY - rect.top;

              var hoveredSpot = null;
              spotData.forEach(function(spot, index) {
                var dist = Math.sqrt((spot.X - x) ** 2 + (y - spot.Y) ** 2);
                if (dist <= 5 && spot.cluster !== 0) {
                  hoveredSpot = spot;
                }
              });

              if (hoveredSpot) {
                $tooltip.text("Domain " + hoveredSpot.cluster);
                $tooltip.show();
                updateTooltipPosition(e);
              } else {
                $tooltip.hide();
              }

              if (drawing && currentCluster > 0) {
                var newSpots = selectSpots(x, y);

                renderCanvas("spotCanvas", true, true);

                newSpots.forEach(function(index) {
                  if (!lastDrawnSpots.has(index)) {
                    var spot = spotData[index];
                    ctx.beginPath();
                    ctx.arc(spot.X, spot.Y, 3, 0, 2 * Math.PI);
                    if (!eraseMode && spot.cluster === 0) {
                      ctx.fillStyle = getClusterColor(currentCluster);
                      ctx.fill();
                    } else if (eraseMode && spot.cluster === currentCluster) {
                      ctx.strokeStyle = "black";
                      ctx.lineWidth = 1;
                      ctx.stroke();
                    }
                    lastDrawnSpots.add(index);
                  }
                });

                ctx.beginPath();
                ctx.arc(x, y, brushSize, 0, 2 * Math.PI);
                ctx.strokeStyle = eraseMode ? "red" : getClusterColor(currentCluster);
                ctx.lineWidth = 1;
                ctx.stroke();
              }
            };

            canvas.onmouseup = function(e) {
              if (currentCluster === 0) return;
              drawing = false;
              lastDrawnSpots = new Set();
              Shiny.setInputValue("selected_spots", Array.from(tempSelectedSpots));
              Shiny.setInputValue("erase_mode", eraseMode);
              renderCanvas("spotCanvas", true, true);
            };
            
            canvas.onmouseenter = function() {
              // 타이머 취소 (마우스가 캔버스에 들어오면 하이라이트 유지)
              if (highlightTimeout) {
                clearTimeout(highlightTimeout);
                highlightTimeout = null;
              }
              renderCanvas("spotCanvas", true, true, null, null);
              Shiny.setInputValue("hovered_cluster", null);
            };
            
            canvas.onmouseleave = function() {
              $tooltip.hide();
            };
            
            function updateTooltipPosition(e) {
              var x = e.pageX + 10;
              var y = e.pageY + 10;
              $tooltip.css({ top: y + "px", left: x + "px" });
            }
          });

          Shiny.addCustomMessageHandler("updateSpotCanvas", function(data) {
            spotData = JSON.parse(JSON.stringify(data.spots));
            currentCluster = data.currentCluster;
            if (data.brushSize) brushSize = data.brushSize;
            clusterInfo = data.clusterInfo || {};
            if (currentCluster !== previousCluster) {
              tempSelectedSpots = new Set();
            }
            previousCluster = currentCluster;
            renderCanvas("spotCanvas", true, true);
          });

          Shiny.addCustomMessageHandler("toggleEraseMode", function(data) {
            eraseMode = data.eraseMode;
            tempSelectedSpots = new Set();
            renderCanvas("spotCanvas", true, true);
          });
          
          Shiny.addCustomMessageHandler("highlightCurrentCluster", function(data) {
            var currentCluster = data.currentCluster;
            // 모든 클러스터 카드에서 current 클래스 제거
            $(".cluster-card-current").removeClass("cluster-card-current");
            // 현재 클러스터 카드에 current 클래스 추가
            $("#cluster-card-" + currentCluster).addClass("cluster-card-current");
          });

          function filterCards(searchTerm, containerId) {
            searchTerm = searchTerm.toLowerCase();
            var cards = document.querySelectorAll(`#${containerId} .gene-card`);
            cards.forEach(function(card) {
              var gene = card.getAttribute("data-gene");
              if (gene.includes(searchTerm)) {
                card.style.display = "";
              } else {
                card.style.display = "none";
              }
            });
          }
        '))
      )
    ),
    tabPanel(
      tags$span(icon("dna"), " Test for DVGs", class = "tab-panel-title"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          step_header("STEP 3. Test for DVGs", "step3"), ## [4] ▶▶
          numericInput("fwer_value", "FWER Threshold:", value = 0.01, min = 0, max = 1, step = "any", width = "100%"),
          tags$small("Any real number between 0 and 1 (exclusive), e.g. 0.01",
                     style = "color: #7f8c8d; display: block; margin-top: -10px; margin-bottom: 12px;"),
          radioButtons("cluster_source", "Tissue Domain Criteria:", 
                       choices = list("Brush-assigned Tissue Domain" = "brush", "Predefined Tissue Domain" = "file"), 
                       selected = "brush", width = "100%"),
          uiOutput("test_button_ui"),
          div(style = "margin-top: 20px;"),
          tags$h4("Plot Adjustment", style = "color: #2980b9; font-size: 16px; margin-bottom: 10px;"),
          sliderInput("test_point_size", "Point Size:", min = 0.1, max = 1, value = 0.2, step = 0.01, width = "100%", ticks = FALSE),
          sliderInput("test_point_alpha", "Point Transparency:", min = 0.1, max = 1, value = 0.8, step = 0.05, width = "100%", ticks = FALSE)
        ),
        mainPanel(
          width = 9,
          div(
            class = "test-results-container",
            style = "padding: 20px; background-color: #fff; border-radius: 8px; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); margin-top: -20px;",
            fluidRow(
              column(
                4,
                tabsetPanel(
                  id = "test_result_tabs",
                  type = "tabs",
                  tabPanel("DVGs",
                           uiOutput("levene_fwer_cards_ui"),
                           div(style = "margin-top: 10px;", uiOutput("download_levene_fwer_ui"))
                  ),
                  tabPanel("DVG-only genes",
                           uiOutput("levene_unique_cards_ui"),
                           div(style = "margin-top: 10px;", uiOutput("download_levene_unique_ui"))
                  )
                )
              ),
              column(
                8,
                uiOutput("test_result_gene_analysis")
              )
            )
          )
        )
      )
    ),
    
    tabPanel(
      tags$span(icon("braille"), " Spatial Domain Variability Testing", class = "tab-panel-title"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          step_header("STEP 4. Spatial Domain Variability Test", "step4"), ## [4] ▶▶
          numericInput("posthoc_alpha", "Significance level (Friedman / Nemenyi):",
                       value = 0.01, min = 0, max = 1, step = "any", width = "100%"),
          tags$small("Any real number between 0 and 1 (exclusive), e.g. 0.01",
                     style = "color: #7f8c8d; display: block; margin-top: -10px; margin-bottom: 12px;"),
          uiOutput("run_tissue_analysis_ui"),
          conditionalPanel(
            condition = "!output.tissue_analysis_ready",
            div(
              style = "color: #7f8c8d; font-style: italic; padding: 10px;",
              "Run STEP 3 (Test for DVGs) first."
            )
          ),
          div(style = "margin-top: 20px;"),
          uiOutput("tissue_friedman_result_sidebar")
        ),
        mainPanel(
          width = 9,
          div(
            class = "tissue-analysis-container",
            style = "padding: 20px; background-color: #fff; border-radius: 8px; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); min-height: 1000px;",
            # Main content row with canvas and variability ranking
            fluidRow(
              column(
                7,
                div(
                  style = "display: flex; flex-direction: column; align-items: flex-start;",
                  div(
                    style = "position: relative; width: 500px; margin-left: 10px;",
                    h3("Tissue Domain Variability", style = "color: #2c3e50; font-size: 20px; position: absolute; top: -30px; left: 50%; transform: translateX(-50%); margin: 0;"),
                    tags$canvas(id = "tissueCanvas", width = "500", height = "500", class = "canvas-container")
                  )
                )
              ),
              column(
                5,
                # Variability Ranking aligned with canvas height
                div(
                  style = "height: 500px; padding-left: 20px; display: flex; flex-direction: column;",
                  # Title at top (aligned with canvas top)
                  div(
                    style = "flex: 0 0 auto; margin-top: -30px; margin-bottom: 10px;",
                    h4("Variability Ranking", style = "color: #2980b9; font-size: 18px; margin: 0;")
                  ),
                  # Cluster cards container - fills remaining space to match canvas height
                  div(
                    style = "flex: 1; max-width: 300px; height: calc(500px - 40px);",
                    uiOutput("tissue_cluster_status_section_cards_only")
                  )
                )
              )
            ),
            # Bottom section with Post-hoc Analysis only
            div(
              style = "margin-top: 20px;",
              fluidRow(
                column(
                  12,
                  # Post-hoc Analysis at bottom
                  conditionalPanel(
                    condition = "output.tissue_analysis_ready",
                    div(
                      h4("Post-hoc Analysis (Nemenyi)", style = "color: #2980b9; margin-bottom: 15px;"),
                      plotOutput("tissue_critical_diff_diagram", height = "180px")
                    )
                  )
                )
              )
            )
          )
        )
      ),
      # Add JavaScript for tissue canvas
      tags$head(
        tags$script(HTML('
      var tissueSpotData = [];
      var tissueClusterRanking = {};
      var tissueHighlightTimeout = null;
      
      function renderTissueCanvas(clear = true, highlightedCluster = null) {
        var canvas = document.getElementById("tissueCanvas");
        if (!canvas) return;
        
        var ctx = canvas.getContext("2d");
        if (clear) {
          ctx.clearRect(0, 0, canvas.width, canvas.height);
        }
        drawTissueSpots(ctx, highlightedCluster);
      }
      
      function drawTissueSpots(ctx, highlightedCluster) {
        tissueSpotData.forEach(function(spot) {
          ctx.beginPath();
          ctx.arc(spot.X, spot.Y, 3, 0, 2 * Math.PI);
          
          var clusterNum = Number(spot.cluster);
          
          if (highlightedCluster !== null) {
            var highlightedNum = Number(highlightedCluster);
            if (clusterNum === highlightedNum && clusterNum !== 0) {
              // 선택된 클러스터는 완전 불투명하게
              ctx.fillStyle = spot.color;
              ctx.globalAlpha = 1.0;
              ctx.fill();
            } else if (clusterNum !== 0) {
              // 다른 클러스터는 반투명하게
              ctx.fillStyle = spot.color;
              ctx.globalAlpha = 0.3;
              ctx.fill();
            } else {
              // 클러스터 0은 테두리만
              ctx.strokeStyle = "black";
              ctx.lineWidth = 1;
              ctx.globalAlpha = 0.3;
              ctx.stroke();
            }
          } else if (clusterNum !== 0) {
            // 하이라이트가 없을 때는 모든 클러스터 동일하게
            ctx.fillStyle = spot.color;
            ctx.globalAlpha = 0.8;
            ctx.fill();
          } else {
            // 클러스터 0은 테두리만
            ctx.strokeStyle = "black";
            ctx.lineWidth = 1;
            ctx.globalAlpha = 0.8;
            ctx.stroke();
          }
          
          // globalAlpha 리셋
          ctx.globalAlpha = 1.0;
        });
      }
      
      Shiny.addCustomMessageHandler("initTissueCanvas", function(data) {
        tissueSpotData = data.spots;
        tissueClusterRanking = data.clusterRanking;
        renderTissueCanvas(true);
        
        var canvas = document.getElementById("tissueCanvas");
        if (!canvas) return;
        
        var $tooltip = $("#customTooltip");
        
        canvas.onmousemove = function(e) {
          var rect = canvas.getBoundingClientRect();
          var x = e.clientX - rect.left;
          var y = e.clientY - rect.top;
          
          var hoveredSpot = null;
          tissueSpotData.forEach(function(spot) {
            var dist = Math.sqrt((spot.X - x) ** 2 + (spot.Y - y) ** 2);
            if (dist <= 5 && spot.cluster !== 0) {
              hoveredSpot = spot;
            }
          });
          
          if (hoveredSpot) {
            var rankInfo = tissueClusterRanking[hoveredSpot.cluster] || {};
            var tooltipText = "Domain " + hoveredSpot.cluster;
            if (rankInfo.rank) {
              tooltipText += " (Rank: " + rankInfo.rank + ")";
            }
            $tooltip.text(tooltipText);
            $tooltip.show();
            var tooltipX = e.pageX + 10;
            var tooltipY = e.pageY + 10;
            $tooltip.css({ top: tooltipY + "px", left: tooltipX + "px" });
          } else {
            $tooltip.hide();
          }
        };
        
        canvas.onmouseleave = function() {
          $tooltip.hide();
          if (tissueHighlightTimeout) {
            clearTimeout(tissueHighlightTimeout);
            tissueHighlightTimeout = null;
          }
          renderTissueCanvas(true, null);
        };
      });
      
      Shiny.addCustomMessageHandler("highlightTissueCluster", function(data) {
        if (tissueHighlightTimeout) {
          clearTimeout(tissueHighlightTimeout);
          tissueHighlightTimeout = null;
        }
        
        renderTissueCanvas(true, data.highlightedCluster);
        
        if (data.highlightedCluster !== null) {
          tissueHighlightTimeout = setTimeout(function() {
            renderTissueCanvas(true, null);
            tissueHighlightTimeout = null;
          }, 2000);
        }
      });
    '))
      )
    ),
    
    tabPanel(
      tags$span(icon("hubspot"), " Spatial Gene Co-expression Analysis", class = "tab-panel-title"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          step_header("STEP 5. Spatial Gene Co-expression Analysis", "step5"), ## [4] ▶▶
          uiOutput("top_n_levene_ui"),
          numericInput("min_cluster_size", "Min Module Size:", value = 5, min = 1, max = 100, step = 1, width = "100%"),
          uiOutput("run_spatial_analysis_ui"),
          conditionalPanel(
            condition = "!output.spatial_analysis_ready",
            div(
              style = "color: #7f8c8d; font-style: italic; padding: 10px;",
              "Run STEP 3 (Test for DVGs) first."
            )
          ),
          div(style = "margin-top: 20px;"),
          uiOutput("gene_search_ui"),
          uiOutput("cluster_table_title"),
          uiOutput("cluster_genes_cards_ui")
        ),
        mainPanel(
          width = 9,
          div(
            class = "spatial-analysis-container",
            style = "padding: 20px; background-color: #fff; border-radius: 8px; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);",
            fluidRow(
              column(
                6,
                visNetworkOutput("cluster_network_plot", height = "500px"),
                div(
                  style = "margin-top: 10px; margin-bottom: 10px;",
                  uiOutput("download_cluster_network_ui")
                ),
                fluidRow(
                  style = "margin-top: 10px;",
                  column(
                    6,
                    sliderInput(
                      "spatial_point_size",
                      "Point Size:",
                      min = 0.1,
                      max = 1,
                      value = 0.2,
                      step = 0.01,
                      width = "100%",
                      ticks = FALSE
                    )
                  ),
                  column(
                    6,
                    sliderInput(
                      "spatial_point_alpha",
                      "Point Transparency:",
                      min = 0.1,
                      max = 1,
                      value = 0.8,
                      step = 0.05,
                      width = "100%",
                      ticks = FALSE
                    )
                  )
                )
              ),
              column(
                6,
                uiOutput("levene_gene_analysis"),
                uiOutput("download_spatial_analysis_ui")
              )
            )
          )
        )
      )
    ), ## [5] ▶▶
    help_tab_ui() ## [5] ▶▶
  )
)


# ============================================================================
# 공통 함수들을 먼저 정의
# ============================================================================

# 클러스터 색상 생성 함수
generate_cluster_colors <- function(clusters) {
  cluster_levels <- levels(clusters)
  non_zero_levels <- cluster_levels[cluster_levels != "0" & cluster_levels != "Non-Cluster"]
  n_colors <- length(non_zero_levels)
  colors <- hcl.colors(n_colors, palette = "Set2")
  setNames(c("#bdc3c7", colors), c("0", non_zero_levels))
}

# 공통 범례 생성 함수
create_cluster_legend <- function(clusters, max_per_row = 5) {
  cluster_levels <- levels(clusters)
  non_zero_levels <- cluster_levels[cluster_levels != "0" & cluster_levels != "Non-Cluster"]
  n_colors <- length(non_zero_levels)
  colors <- hcl.colors(n_colors, palette = "Set2")
  
  items <- lapply(1:n_colors, function(i) {
    tags$span(
      style = "display: inline-block; margin-right: 10px;",
      tags$span(
        style = sprintf("display: inline-block; width: 12px; height: 12px; background-color: %s; border-radius: 50%%; margin-right: 5px;", colors[i])
      ),
      tags$span(style = "font-size: 12px;", non_zero_levels[i])
    )
  })
  
  rows <- split(items, ceiling(seq_along(items) / max_per_row))
  div(
    style = "text-align: left; margin: 0; padding: 0;",
    lapply(rows, function(row) div(style = "margin-bottom: 4px;", row))
  )
}

# 공통 발현 범례 생성 함수
create_expression_legend <- function(expr_range) {
  min_val <- expr_range[1]
  max_val <- expr_range[2]
  gradient_style <- "background: linear-gradient(to right, #440154FF, #414487FF, #2A8A8CFF, #5EC962FF, #FDE725FF);"
  
  tagList(
    div(
      style = "text-align: right; margin: 0; padding: 0;",
      div(style = "font-size: 12px; margin-bottom: 4px;", "Expression Level"),
      div(
        style = "position: relative; width: 150px; margin-left: auto;",
        div(style = paste0(gradient_style, "width: 150px; height: 12px; border-radius: 2px;")),
        div(
          style = "font-size: 12px; margin-top: 4px; position: relative; width: 150px;",
          span(sprintf("%.2f", min_val), style = "position: absolute; left: 0; color: #2c3e50;"),
          span(sprintf("%.2f", max_val), style = "position: absolute; right: 0; color: #2c3e50;")
        )
      )
    )
  )
}

# 공통 클러스터 분포 플롯 함수
create_cluster_plot <- function(coords, clusters, point_size, point_alpha, title = "Tissue Domains") {
  levels(clusters) <- ifelse(levels(clusters) == "0", "Non-Cluster", levels(clusters))
  df <- data.frame(x = coords$X, y = coords$Y, Cluster = clusters)
  color_mapping <- generate_cluster_colors(clusters)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data = df[df$Cluster == "Non-Cluster", ], 
               color = "#bdc3c7", shape = 1, size = 5 * point_size, alpha = point_alpha) +
    geom_point(data = df[df$Cluster != "Non-Cluster", ], 
               aes(color = Cluster), shape = 16, size = 5 * point_size, alpha = point_alpha) +
    scale_color_manual(values = color_mapping) +
    labs(title = title, x = "X Coordinate", y = "Y Coordinate") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

# 공통 발현 플롯 함수
create_expression_plot <- function(coords, clusters, expression, gene, point_size, point_alpha) {
  levels(clusters) <- ifelse(levels(clusters) == "0", "Non-Cluster", levels(clusters))
  df <- data.frame(x = coords$X, y = coords$Y, Expression = expression, Cluster = clusters)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data = df[df$Cluster == "Non-Cluster", ], 
               color = "#bdc3c7", shape = 1, size = 5 * point_size, alpha = point_alpha) +
    geom_point(data = df[df$Cluster != "Non-Cluster", ], 
               aes(color = Expression), shape = 16, size = 5 * point_size, alpha = point_alpha) +
    scale_color_viridis_c(option = "viridis", name = "Expression") +
    labs(title = paste("Expression of", gene), x = "X Coordinate", y = "Y Coordinate") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

# 과학적 표기법 형식화 함수
format_scientific_notation <- function(value) {
  if (length(value) == 0 || is.na(value)) return("NA")
  sci <- format(value, digits = 5, scientific = TRUE)
  parts <- strsplit(sci, "e")[[1]]
  mantissa <- parts[1]
  exponent <- as.numeric(parts[2])
  superscript_digits <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
  superscript_minus <- "⁻"
  
  if (is.na(exponent)) return(paste0(mantissa, " x 10^NA"))
  
  exp_digits <- strsplit(as.character(abs(exponent)), "")[[1]]
  exp_str <- if (exponent < 0) {
    paste0(superscript_minus, paste(sapply(exp_digits, function(d) {
      idx <- as.numeric(d) + 1
      if (idx >= 1 && idx <= length(superscript_digits)) superscript_digits[idx] else d
    }), collapse = ""))
  } else {
    paste(sapply(exp_digits, function(d) {
      idx <- as.numeric(d) + 1
      if (idx >= 1 && idx <= length(superscript_digits)) superscript_digits[idx] else d
    }), collapse = "")
  }
  
  paste0(mantissa, " x 10", exp_str)
}

## [3] ▶
scale_spots_to_canvas <- function(coords, canvas_width = 500, canvas_height = 500) {
  spots <- coords
  spots$cluster <- 0
  max_x <- max(spots$X, na.rm = TRUE)
  max_y <- max(spots$Y, na.rm = TRUE)
  scale_factor <- min(canvas_width / max_x, canvas_height / max_y) * 0.9
  spots$X <- spots$X * scale_factor + (canvas_width  - max_x * scale_factor) / 2
  spots$Y <- spots$Y * scale_factor + (canvas_height - max_y * scale_factor) / 2
  spots$Y <- canvas_height - spots$Y   # Y축 반전 (기존 코드와 동일)
  spots
}

spots_to_list <- function(spots) {
  has_uploaded <- "uploadedCluster" %in% colnames(spots)
  lapply(seq_len(nrow(spots)), function(i) {
    item <- list(X = spots$X[i], Y = spots$Y[i], cluster = spots$cluster[i])
    if (has_uploaded) item$uploadedCluster <- spots$uploadedCluster[i]
    item
  })
}

# 0 < x < 1 인 실수인지 확인. 아니면 에러 알림을 띄우고 FALSE 반환
check_open_unit_interval <- function(x, label) {
  if (is.null(x) || is.na(x) || x <= 0 || x >= 1) {
    showNotification(
      sprintf("%s must be a number strictly between 0 and 1 (e.g., 0.01).", label),
      type = "error", duration = 8
    )
    return(FALSE)
  }
  TRUE
}

# 서버 정의
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 1024*1024^2)
  
  ## [6] ▶▶
  ### Help: STEP 제목 옆 '?' -> 모달, 모달의 'Open full tutorial' -> Help 탭으로 이동
  for (s in c("step1", "step2", "step3", "step4", "step5")) {
    local({
      sec <- s
      observeEvent(input[[paste0("help_", sec)]], help_modal(sec))
    })
  }
  observeEvent(input$goto_help_tab, {
    removeModal()
    updateTabsetPanel(session, "main_tabs", selected = "help_tab")
  })
  
  expr_data <- reactiveVal(NULL)
  coords_data <- reactiveVal(NULL)
  spots_data <- reactiveVal(NULL)
  brush_clusters <- reactiveVal(NULL)
  current_cluster <- reactiveVal(0)
  fwer_results <- reactiveVal(NULL)
  cluster_info <- reactiveVal(list())
  uploaded_clusters <- reactiveVal(NULL)
  erase_mode <- reactiveVal(FALSE)
  data_matched <- reactiveVal(TRUE)
  
  selected_gene <- reactiveVal(NULL)
  selected_gene_stats <- reactiveVal(NULL)
  
  # 데이터 크기 일치 여부 확인, 툴팁 메시지 설정, 초기화
  observe({
    expr <- expr_data()
    coords <- coords_data()
    files_present <- !is.null(expr) && !is.null(coords)
    
    # 데이터 상태 확인 및 초기화
    if (!files_present || (files_present && ncol(expr) != nrow(coords))) {
      # 데이터가 하나라도 없거나 크기가 일치하지 않으면 초기화
      data_matched(FALSE)
      fwer_results(NULL)           # Test Results 초기화
      selected_gene(NULL)         # 클릭한 유전자 초기화
      selected_gene_stats(NULL)   # 유전자 통계 초기화
      spots_data(NULL)            # Canvas 데이터 초기화
      brush_clusters(NULL)        # Canvas 클러스터 초기화
      current_cluster(0)          # 현재 클러스터 초기화
      cluster_info(list())        # 클러스터 정보 초기화
      
      # Spatial Analysis 초기화
      levene_graph_data(NULL)         # 그래프 데이터 초기화
      selected_cluster_node(NULL)     # 선택된 클러스터 노드 초기화
      selected_levene_gene(NULL)      # 선택된 Levene 유전자 초기화
      selected_levene_gene_stats(NULL) # Levene 유전자 통계 초기화
      
      # 캔버스 초기화 메시지 전송
      session$sendCustomMessage("initDataCanvas", list(
        spots = list(),
        uploadedClusters = NULL
      ))
      session$sendCustomMessage("initSpotCanvas", list(
        spots = list(),
        currentCluster = 0,
        brushSize = input$brush_size,
        clusterInfo = list(),
        uploadedClusters = NULL
      ))
      
      # 모든 버튼에 대해 툴팁 메시지 설정
      buttons <- c("test_button", "new_cluster_btn", "run_spatial_analysis")
      for (button_id in buttons) {
        if (!files_present) {
          session$sendCustomMessage("updateTooltipMessage", list(buttonId = button_id, text = "Please load your data"))
        } else {
          session$sendCustomMessage("updateTooltipMessage", list(buttonId = button_id, text = "Expression data doesn't match with coordinate data"))
        }
      }
      
      if (files_present && ncol(expr) != nrow(coords)) {
        showNotification(
          paste("Number of spots (", nrow(coords), ") does not match number of samples in expression matrix (", ncol(expr), ")."),
          type = "error",
          duration = 10
        )
      }
    } else {
      # 데이터가 모두 있고 크기가 일치하면 data_matched TRUE로 설정
      data_matched(TRUE)
      
      # 모든 버튼에 대해 툴팁 메시지 비우기
      buttons <- c("test_button", "gene_button", "new_cluster_btn", "run_spatial_analysis")
      for (button_id in buttons) {
        session$sendCustomMessage("updateTooltipMessage", list(buttonId = button_id, text = ""))
      }
    }
  })
  
  ### STEP 1: 데이터 업로드 - 좌표 데이터 처리
  observeEvent(input$coords_file, {
    req(input$coords_file)
    withProgress(message = "Loading Coordinates...", value = 0, {
      incProgress(0.5, detail = "Loading coordinates...")
      coords_ext <- tools::file_ext(input$coords_file$datapath)
      if (coords_ext != "csv") {
        showNotification("Unsupported coordinates file type. Please use .csv.", type = "error")
        return()
      }
      
      coords <- read.csv(input$coords_file$datapath, check.names = FALSE)
      
      # 열 개수 확인
      if (ncol(coords) < 3) {
        showNotification("Coordinates file must have at least 3 columns (Spot, X, Y).", type = "error")
        return()
      }
      coords <- coords[, 1:3]  # 첫 3열만 사용 (Spot, X, Y)
      colnames(coords) <- c("Spot", "X", "Y")
      
      # NA 또는 비수치형 값 확인 및 필터링
      if (any(is.na(coords$X)) || any(is.na(coords$Y))) {
        showNotification("Warning: Some X or Y values are NA or non-numeric and will be filtered out.", type = "warning")
      }
      coords <- coords %>% filter(!is.na(X) & !is.na(Y))
      if (nrow(coords) == 0) {
        showNotification("No valid numeric X, Y coordinates found after filtering.", type = "error")
        return()
      }
      
      # 최소값 빼기 (원점으로 바꾸기 위해)
      coords$X <- coords$X - min(coords$X, na.rm = TRUE)
      coords$Y <- coords$Y - min(coords$Y, na.rm = TRUE)
      coords_data(coords)
      
      # 데이터 초기화
      spots <- coords
      spots$cluster <- 0
      brush_clusters(rep(0, nrow(spots)))
      spots_data(spots)
      current_cluster(0)
      cluster_info(list())
      erase_mode(FALSE)
      uploaded_clusters(NULL)  # Cluster Vector 초기화
      
      expr <- expr_data()
      if (!is.null(expr)) {
        expr_colnames <- colnames(expr)
        if (nrow(spots) != length(expr_colnames)) {
          data_matched(FALSE)
          showNotification(
            paste("Number of spots (", nrow(spots), ") does not match number of samples in expression matrix (", length(expr_colnames), ")."),
            type = "error",
            duration = 10
          )
          return()  # 크기 불일치 시 여기서 중단
        } else {
          data_matched(TRUE)
          rownames(spots) <- expr_colnames
        }
      }
      
      # uploaded_clusters가 유효할 경우에만 추가
      uploaded <- uploaded_clusters()
      if (!is.null(uploaded) && length(uploaded) == nrow(spots)) {
        spots$uploadedCluster <- uploaded
      }
      
      canvas_width <- 500
      canvas_height <- 500
      max_x <- max(spots$X, na.rm = TRUE)
      max_y <- max(spots$Y, na.rm = TRUE)
      scale_factor <- min(canvas_width / max_x, canvas_height / max_y) * 0.9
      spots$X <- spots$X * scale_factor
      spots$Y <- spots$Y * scale_factor
      offset_x <- (canvas_width - max_x * scale_factor) / 2
      offset_y <- (canvas_height - max_y * scale_factor) / 2
      spots$X <- spots$X + offset_x
      spots$Y <- spots$Y + offset_y
      spots$Y <- canvas_height - spots$Y  # Y축 반전을 offset 적용 후에
      
      spots_data(spots)
      
      incProgress(1.0, detail = "Rendering canvas...")
      spots_list <- spots_to_list(spots)
      
      # uploadedClusters가 있을 때만 전달
      uploaded_clusters_data <- if ("uploadedCluster" %in% colnames(spots)) spots$uploadedCluster else NULL
      
      session$sendCustomMessage("initDataCanvas", list(
        spots = spots_list,
        uploadedClusters = uploaded_clusters_data
      ))
      
      session$sendCustomMessage("initSpotCanvas", list(
        spots = spots_list,
        currentCluster = isolate(current_cluster()),
        brushSize = input$brush_size,
        clusterInfo = cluster_info(),
        uploadedClusters = uploaded_clusters_data
      ))
    })
  })
  
  ### Expression Matrix 처리
  observeEvent(input$expr_file, {
    req(input$expr_file)
    withProgress(message = "Loading Expression Matrix...", value = 0, {
      incProgress(0.5, detail = "")
      expr_ext <- tools::file_ext(input$expr_file$datapath)
      if (expr_ext != "csv") {
        showNotification("Unsupported expression file type. Please use .csv.", type = "error")
        return()
      }
      
      expr <- read.csv(input$expr_file$datapath, row.names = 1, check.names = FALSE)
      
      # 모든 데이터 초기화
      expr_data(as.matrix(expr))
      spots_data(NULL)
      brush_clusters(NULL)
      current_cluster(0)
      cluster_info(list())
      erase_mode(FALSE)
      uploaded_clusters(NULL)  # Cluster Vector 초기화
      
      coords <- coords_data()
      if (!is.null(coords)) {
        spots <- coords
        spots$cluster <- 0
        brush_clusters(rep(0, nrow(spots)))
        expr_colnames <- colnames(expr)
        if (nrow(spots) != length(expr_colnames)) {
          data_matched(FALSE)
          showNotification(
            paste("Number of spots (", nrow(spots), ") does not match number of samples in expression matrix (", length(expr_colnames), ")."),
            type = "error",
            duration = 10
          )
          return()  # 크기 불일치 시 여기서 중단
        } else {
          data_matched(TRUE)
          rownames(spots) <- expr_colnames
        }
        
        canvas_width <- 500
        canvas_height <- 500
        max_x <- max(spots$X, na.rm = TRUE)
        max_y <- max(spots$Y, na.rm = TRUE)
        scale_factor <- min(canvas_width / max_x, canvas_height / max_y) * 0.9
        spots$X <- spots$X * scale_factor
        spots$Y <- spots$Y * scale_factor
        offset_x <- (canvas_width - max_x * scale_factor) / 2
        offset_y <- (canvas_height - max_y * scale_factor) / 2
        spots$X <- spots$X + offset_x
        spots$Y <- spots$Y + offset_y
        spots$Y <- canvas_height - spots$Y  # Y축 반전을 offset 적용 후에
        spots_data(spots)
        
        spots_list <- spots_to_list(spots)
        
        uploaded_clusters_data <- if ("uploadedCluster" %in% colnames(spots)) spots$uploadedCluster else NULL
        
        session$sendCustomMessage("initDataCanvas", list(
          spots = spots_list,
          uploadedClusters = uploaded_clusters_data
        ))
        
        session$sendCustomMessage("initSpotCanvas", list(
          spots = spots_list,
          currentCluster = isolate(current_cluster()),
          brushSize = input$brush_size,
          clusterInfo = cluster_info(),
          uploadedClusters = uploaded_clusters_data
        ))
      }
      
      incProgress(1.0, detail = "Expression matrix loaded")
    })
  })
  
  output$cluster_file_ui <- renderUI({
    is_data_matched <- data_matched()
    expr <- expr_data()
    coords <- coords_data()
    files_present <- !is.null(expr) && !is.null(coords)
    
    input_field <- fileInput(
      "cluster_file",
      HTML("Predefined Tissue Domain<br>(Optional; .csv)"),
      accept = ".csv",
      width = "100%"
    )
    
    if (!files_present || !is_data_matched) {
      div(
        input_field,
        style = "pointer-events: none; opacity: 0.5;"
      )
    } else {
      input_field
    }
  })
  
  ### 클러스터 파일 업로드 처리
  observeEvent(input$cluster_file, {
    req(input$cluster_file)
    withProgress(message = "Loading Domain Labels...", value = 0, {
      incProgress(0.5, detail = "Processing domain label file...")
      cluster_ext <- tools::file_ext(input$cluster_file$datapath)
      
      if (cluster_ext != "csv") {
        showNotification("Unsupported domain label file type. Please use .csv.", type = "error")
        return()
      }
      
      clusters_df <- read.csv(input$cluster_file$datapath, header = TRUE, check.names = FALSE)
      
      # Check if clusters_df is a data frame and has data
      if (is.null(clusters_df) || !is.data.frame(clusters_df) || nrow(clusters_df) == 0) {
        showNotification("Domain label file is empty or could not be read properly.", type = "error")
        return()
      }
      
      # Check if it has exactly 2 columns
      if (ncol(clusters_df) != 2) {
        showNotification("Domain label file must contain exactly 2 columns: spot ID and domain label.", type = "error")
        return()
      }
      
      # Assign column names if not already present
      if (!all(c("Spot", "Cluster") %in% colnames(clusters_df))) {
        colnames(clusters_df) <- c("Spot", "Cluster")
      }
      
      # Extract spot names and cluster values
      spot_names <- clusters_df[["Spot"]]
      clusters <- clusters_df[["Cluster"]]
      
      # Validate spot names and clusters
      if (is.null(spot_names) || length(spot_names) == 0) {
        showNotification("No valid spot names found in the first column.", type = "error")
        return()
      }
      
      if (is.null(clusters) || length(clusters) == 0) {
        showNotification("No valid domain labels found in the second column.", type = "error")
        return()
      }
      
      uploaded_clusters(clusters)
      
      # Process spots_data if available
      spots <- spots_data()
      if (is.null(spots)) {
        showNotification("Spot data is not loaded yet. Please load coordinates first.", type = "warning")
      } else {
        if (length(clusters) != nrow(spots) || length(spot_names) != nrow(spots)) {
          showNotification(
            paste("Domain label file rows (", length(clusters), ") do not match the number of spots (", nrow(spots), ")."),
            type = "error",
            duration = 10
          )
          return()
        }
        
        if (!all(spot_names %in% spots$Spot)) {
          showNotification("Spot IDs in the domain label file do not match those in the coordinate file.", type = "warning")
        }
        
        spots$uploadedCluster <- clusters
        spots_data(spots)
        
        spots_list <- spots_to_list(spots)
        
        session$sendCustomMessage("initDataCanvas", list(
          spots = spots_list,
          uploadedClusters = spots$uploadedCluster
        ))
        
        session$sendCustomMessage("initSpotCanvas", list(
          spots = spots_list,
          currentCluster = isolate(current_cluster()),
          brushSize = input$brush_size,
          clusterInfo = cluster_info()
        ))
      }
      
      incProgress(1.0, detail = "Domain labels loaded")
    })
  })
  
  ## [4] ▶
  ### Example data 로드 (P16 melanoma brain metastasis + BANYAN domains)
  observeEvent(input$load_example_btn, {
    if (is.null(example_data)) {
      showNotification("Example data files were not found on the server (test_data/ folder).",
                       type = "error", duration = 10)
      return()
    }
    
    withProgress(message = "Loading Example Data...", value = 0, {
      incProgress(0.2, detail = example_data$description)
      
      expr   <- example_data$expr
      coords <- example_data$coords
      labels <- example_data$labels
      
      ## 1) 좌표 원점 이동 (coords_file observer 와 동일)
      coords$X <- coords$X - min(coords$X, na.rm = TRUE)
      coords$Y <- coords$Y - min(coords$Y, na.rm = TRUE)
      
      ## 2) 이전 세션 결과 초기화 (파일 업로드 경로와 동일하게 downstream 전부 리셋)
      fwer_results(NULL); selected_gene(NULL); selected_gene_stats(NULL)
      tissue_analysis_results(NULL); tissue_cluster_data(NULL)
      levene_graph_data(NULL); selected_cluster_node(NULL)
      selected_levene_gene(NULL); selected_levene_gene_stats(NULL)
      current_cluster(0); cluster_info(list()); erase_mode(FALSE)
      
      ## 3) 데이터 세팅 — 파일 3개를 업로드한 것과 같은 상태
      expr_data(expr)
      coords_data(coords)
      uploaded_clusters(labels)
      data_matched(TRUE)
      
      incProgress(0.5, detail = "Rendering canvas...")
      spots <- scale_spots_to_canvas(coords)
      rownames(spots) <- colnames(expr)
      spots$uploadedCluster <- labels
      brush_clusters(rep(0, nrow(spots)))
      spots_data(spots)
      
      spots_list <- spots_to_list(spots)
      session$sendCustomMessage("initDataCanvas", list(
        spots = spots_list,
        uploadedClusters = labels
      ))
      session$sendCustomMessage("initSpotCanvas", list(
        spots = spots_list,
        currentCluster = 0,
        brushSize = input$brush_size,
        clusterInfo = list(),
        uploadedClusters = labels
      ))
      session$sendCustomMessage("toggleEraseMode", list(eraseMode = FALSE))
      
      ## 4) BANYAN 라벨이 들어왔으므로 DVG 탭의 'Tissue Domain Criteria' 를 Predefined 로 전환
      updateRadioButtons(session, "cluster_source", selected = "file")
      
      incProgress(1.0, detail = "Done")
    })
    
    showNotification(
      paste0("Example data loaded (", example_data$description, "). ",
             "BANYAN tissue domains are preselected \u2014 open 'Test for DVGs' and press 'Test Start!'."),
      type = "message", duration = 8
    )
  })
  
  
  ### STEP 2: 클러스터링 (브러쉬 클러스터만 사용)
  
  output$new_cluster_btn_ui <- renderUI({
    current_info <- cluster_info()
    current_cluster_key <- as.character(current_cluster())
    has_spots <- current_cluster() == 0 || (!is.null(current_info[[current_cluster_key]]) && length(current_info[[current_cluster_key]]) > 0)
    
    is_data_loaded <- !is.null(expr_data()) && !is.null(coords_data())
    is_data_matched <- data_matched()
    files_present <- !is.null(expr_data()) && !is.null(coords_data())
    
    # 버튼 활성화 조건
    is_enabled <- is_data_loaded && is_data_matched && (current_cluster() == 0 || has_spots)
    
    # 버튼 정의
    if (current_cluster() == 0) {
      btn <- actionButton(
        "new_cluster_btn",
        tags$span(icon("play"), "Start!"),
        class = "btn-primary btn-sm",
        width = "100%"
      )
    } else {
      btn <- actionButton(
        "new_cluster_btn",
        tags$span(icon("plus"), " New Domain"),
        class = "btn-primary btn-sm",
        width = "100%"
      )
    }
    
    if (!files_present || !is_data_matched) {
      btn <- tagAppendAttributes(btn, disabled = "disabled")
    }
    
    btn
  })
  
  # 동적 Erase/Brush 버튼 렌더링
  output$erase_btn_ui <- renderUI({
    if (erase_mode()) {
      actionButton(
        "erase_btn",
        tags$span(icon("paint-brush"), " Brush"),
        class = "btn-primary btn-sm",
        width = "100%"
      )
    } else {
      actionButton(
        "erase_btn",
        tags$span(icon("eraser"), " Erase"),
        class = "btn-warning btn-sm",
        width = "100%"
      )
    }
  })
  
  # Cluster Status 섹션
  output$cluster_status_section <- renderUI({
    if (is.null(expr_data()) || is.null(coords_data())) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; padding: 10px; text-align: center; white-space: nowrap;",
        "Please upload data to view domain status"
      ))
    }
    
    tagList(
      div(textOutput("current_cluster"), style = "font-size: 14px; margin-top: 10px; margin-bottom: 5px; color: #6c757d;"),
      h4("Domain Status", style = "color: #2980b9; font-size: 18px; margin-top: 0; margin-bottom: 10px;"),
      uiOutput("cluster_status_cards_ui")
    )
  })
  
  output$cluster_status_cards_ui <- renderUI({
    spots <- spots_data()
    if (is.null(spots)) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; text-align: center; padding: 20px;",
        "No domain data available."
      ))
    }
    
    cluster_counts <- table(spots$cluster)
    cluster_df <- data.frame(
      Cluster = names(cluster_counts),
      Spots = as.integer(cluster_counts)
    )
    
    # 클러스터 색상 함수
    getClusterColor <- function(clusterIndex) {
      if (clusterIndex == "0") return("#000000")
      hue <- (as.numeric(clusterIndex) * 137.5) %% 360
      return(paste0("hsl(", hue, ", 70%, 50%)"))
    }
    
    # 카드 리스트 생성
    cards <- lapply(1:nrow(cluster_df), function(i) {
      cluster <- cluster_df$Cluster[i]
      spots <- cluster_df$Spots[i]
      cluster_name <- ifelse(cluster == "0", "Unassigned (0)", paste("Domain", cluster))
      cluster_color <- getClusterColor(cluster)
      
      # 현재 클러스터 여부 확인
      is_current <- as.numeric(cluster) == current_cluster()
      card_class <- if (is_current) "gene-card cluster-card-current" else "gene-card"
      
      div(
        class = card_class,
        id = paste0("cluster-card-", cluster),
        onclick = sprintf("
        var timestamp = new Date().getTime(); 
        Shiny.setInputValue('clicked_cluster_card', '%s_' + timestamp, {priority: 'event'});
      ", cluster),
        style = paste0("border-left: 4px solid ", cluster_color, ";"),
        div(
          class = "gene-card-header", 
          style = paste0("color: ", cluster_color, "; display: flex; align-items: center;"),
          tags$span(
            style = paste0("display: inline-block; width: 12px; height: 12px; border-radius: 50%; background-color: ", cluster_color, "; margin-right: 8px;")
          ),
          cluster_name
        ),
        div(
          class = "gene-card-content",
          tags$span(sprintf("Spots: %d", spots))
        )
      )
    })
    
    # 검색 입력란 없이 카드만 반환
    div(class = "cards-container", id = "cluster-status", cards)
  })
  
  observeEvent(input$clicked_cluster_card, {
    req(spots_data(), input$clicked_cluster_card)
    
    # 타임스탬프가 포함된 클러스터 ID에서 실제 클러스터 번호만 추출
    cluster_id_raw <- input$clicked_cluster_card
    cluster_id <- tryCatch({
      # 언더스코어 이전 부분만 추출 (타임스탬프 제거)
      cluster_part <- strsplit(as.character(cluster_id_raw), "_")[[1]][1]
      as.numeric(cluster_part)
    }, error = function(e) {
      showNotification("Invalid domain ID", type = "error")
      return(NULL)
    })
    
    if (is.null(cluster_id) || is.na(cluster_id)) {
      return()
    }
    
    spots <- spots_data()
    
    spots_list <- spots_to_list(spots)
    
    session$sendCustomMessage("highlightCluster", list(
      spots = spots_list,
      highlightedCluster = cluster_id,  # 숫자로 전달
      currentCluster = current_cluster(),
      brushSize = input$brush_size,
      clusterInfo = cluster_info()
    ))
  })
  
  # Current Cluster 텍스트 출력
  output$current_cluster <- renderText({
    paste("Current Tissue Domain:", current_cluster(), ifelse(erase_mode(), " (Erase Mode)", ""))
  })
  
  observeEvent(input$new_cluster_btn, {
    if (current_cluster() > 0) {
      current_info <- cluster_info()
      current_cluster_key <- as.character(current_cluster())
      if (is.null(current_info[[current_cluster_key]]) || length(current_info[[current_cluster_key]]) == 0) {
        showNotification(
          "At least one spot is needed to create a new domain.",
          type = "warning",
          duration = 5
        )
        return()
      }
    }
    
    current_cluster(current_cluster() + 1)
    erase_mode(FALSE)
    spots <- spots_data()
    
    spots_list <- spots_to_list(spots)
    
    session$sendCustomMessage("updateSpotCanvas", list(
      spots = spots_list,
      currentCluster = current_cluster(),
      clusterInfo = cluster_info()
    ))
    session$sendCustomMessage("toggleEraseMode", list(eraseMode = FALSE))
  })
  
  observeEvent(input$erase_btn, {
    erase_mode(!erase_mode())
    session$sendCustomMessage("toggleEraseMode", list(eraseMode = erase_mode()))
  })
  
  observeEvent(input$selected_spots, {
    spots <- spots_data()
    selected_indices <- input$selected_spots
    if (length(selected_indices) > 0) {
      if (!input$erase_mode) {
        valid_indices <- selected_indices[sapply(selected_indices, function(i) spots$cluster[i + 1] == 0)]
        if (length(valid_indices) > 0) {
          spots$cluster[valid_indices + 1] <- current_cluster()
          brush_clusters(spots$cluster)
          spots_data(spots)
          
          current_info <- cluster_info()
          if (is.null(current_info[[as.character(current_cluster())]])) {
            current_info[[as.character(current_cluster())]] <- valid_indices + 1
          } else {
            current_info[[as.character(current_cluster())]] <- unique(c(current_info[[as.character(current_cluster())]], valid_indices + 1))
          }
          cluster_info(current_info)
        }
      } else {
        valid_indices <- selected_indices[sapply(selected_indices, function(i) spots$cluster[i + 1] == current_cluster())]
        if (length(valid_indices) > 0) {
          spots$cluster[valid_indices + 1] <- 0
          brush_clusters(spots$cluster)
          spots_data(spots)
          
          current_info <- cluster_info()
          current_info[[as.character(current_cluster())]] <- setdiff(current_info[[as.character(current_cluster())]], valid_indices + 1)
          if (length(current_info[[as.character(current_cluster())]]) == 0) {
            current_info[[as.character(current_cluster())]] <- NULL
          }
          cluster_info(current_info)
        }
      }
      
      spots_list <- spots_to_list(spots)
      
      session$sendCustomMessage("updateSpotCanvas", list(
        spots = spots_list,
        currentCluster = current_cluster(),
        brushSize = input$brush_size,
        clusterInfo = cluster_info()
      ))
    }
  })
  
  observeEvent(input$auto_cluster_btn, {
    req(spots_data())
    spots <- spots_data()
    current_info <- cluster_info()
    current_cluster_key <- as.character(current_cluster())
    
    # 남은 스팟(클러스터 0인 스팟) 확인
    remaining <- which(spots$cluster == 0)
    if (length(remaining) == 0) {
      showNotification("No unassigned spots remain.", type = "warning")
      return()
    }
    
    # 현재 클러스터에 스팟이 있는지 확인
    has_spots <- !is.null(current_info[[current_cluster_key]]) && length(current_info[[current_cluster_key]]) > 0
    
    if (has_spots) {
      # 현재 클러스터에 스팟이 있으면 새로운 클러스터 생성
      current_cluster(current_cluster() + 1)
    } else if (current_cluster() == 0) {
      # 클러스터링이 시작되지 않았으면 첫 클러스터 생성
      current_cluster(1)
    }
    # 현재 클러스터가 0이 아니고 스팟이 없으면 현재 클러스터를 유지하고 진행
    
    # 남은 스팟을 현재 클러스터에 할당
    spots$cluster[remaining] <- current_cluster()
    brush_clusters(spots$cluster)
    spots_data(spots)
    
    # 클러스터 정보 업데이트
    current_info[[as.character(current_cluster())]] <- remaining
    cluster_info(current_info)
    
    spots_list <- spots_to_list(spots)
    
    # 캔버스 업데이트
    session$sendCustomMessage("updateSpotCanvas", list(
      spots = spots_list,
      currentCluster = current_cluster(),
      clusterInfo = cluster_info()
    ))
    
    # Erase 모드 비활성화
    erase_mode(FALSE)
    session$sendCustomMessage("toggleEraseMode", list(eraseMode = FALSE))
  })
  
  observeEvent(input$brush_size, {
    req(spots_data())
    spots <- spots_data()
    spots_list <- spots_to_list(spots)
    session$sendCustomMessage("updateSpotCanvas", list(
      spots = spots_list,
      currentCluster = current_cluster(),
      brushSize = input$brush_size,
      clusterInfo = cluster_info()
    ))
  })
  
  ### Reset 버튼 처리
  observeEvent(input$reset_cluster_btn, {
    req(spots_data())
    spots <- spots_data()
    spots$cluster <- 0
    brush_clusters(rep(0, nrow(spots)))
    current_cluster(0)
    cluster_info(list())
    spots_data(spots)
    
    spots_list <- spots_to_list(spots)
    
    session$sendCustomMessage("updateSpotCanvas", list(
      spots = spots_list,
      currentCluster = current_cluster(),
      clusterInfo = cluster_info()
    ))
    erase_mode(FALSE)
    session$sendCustomMessage("toggleEraseMode", list(eraseMode = FALSE))
  })
  
  # Cluster Assignment 또는 Cluster Source 변경 시 Test Results 및 Gene Analysis 초기화
  observeEvent(list(brush_clusters(), input$cluster_source), {
    # Test Results 관련 데이터 초기화
    fwer_results(NULL)
    selected_gene(NULL)
    selected_gene_stats(NULL)
    
  }, ignoreInit = TRUE)
  
  ### STEP 3: 통계 테스트
  observeEvent(input$test_button, {
    if (!check_open_unit_interval(input$fwer_value, "FWER Threshold")) return()
    withProgress(message = "Running Tests...", value = 0, {
      req(expr_data(), spots_data())
      
      incProgress(0.05, detail = "Preparing data...")
      
      expr <- expr_data()
      
      # 클러스터 데이터 준비
      if (input$cluster_source == "brush") {
        clusters_raw <- brush_clusters()
      } else if (input$cluster_source == "file") {
        req(uploaded_clusters())
        clusters_raw <- uploaded_clusters()
        if (length(clusters_raw) != ncol(expr)) {
          stop("Uploaded domain label vector length does not match expression matrix columns.")
        }
      }
      
      # 인덱스 기반 필터링
      selected_indices <- which(clusters_raw != 0)
      
      if (length(selected_indices) < 10) {
        showNotification("Error: Too few spots selected.", type = "error")
        return()
      }
      
      # 행렬 서브셋팅
      expr_selected <- expr[, selected_indices, drop = FALSE]
      clusters <- factor(clusters_raw[selected_indices])
      
      if (nlevels(clusters) < 2) {
        showNotification("Error: At least 2 domains (excluding unassigned spots) are required for test.", type = "error")
        return()
      }
      
      n_genes <- nrow(expr_selected)
      gene_names <- rownames(expr_selected)
      
      incProgress(0.05, detail = sprintf("Testing %d genes (parallel)...", n_genes))
      
      # 병렬 처리 설정
      num_cores <- max(1, detectCores() - 1)
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      # 한 번에 전체 병렬 처리
      results <- foreach(i = 1:n_genes, .combine = 'rbind', .packages = c('car')) %dopar% {
        gene_data <- expr_selected[i, ]
        gene_var <- var(gene_data, na.rm = TRUE)
        
        if (is.na(gene_var) || gene_var == 0) {
          return(data.frame(Gene = gene_names[i], anova_p = NA, levene_p = NA))
        }
        
        anova_p <- tryCatch({
          aov_res <- aov(gene_data ~ clusters)
          summary(aov_res)[[1]]$"Pr(>F)"[1]
        }, error = function(e) NA)
        
        levene_p <- tryCatch({
          lev_res <- leveneTest(gene_data, clusters, center = median)
          lev_res$"Pr(>F)"[1]
        }, error = function(e) NA)
        
        data.frame(Gene = gene_names[i], anova_p = anova_p, levene_p = levene_p)
      }
      
      stopCluster(cl)
      
      incProgress(0.70, detail = "Calculating adjusted p-values...")
      
      # 결과 정리
      anova_df <- data.frame(Gene = results$Gene, P_value = results$anova_p)
      levene_df <- data.frame(Gene = results$Gene, P_value = results$levene_p)
      anova_df$`Adjusted P-value` <- p.adjust(anova_df$P_value, method = "bonferroni")
      levene_df$`Adjusted P-value` <- p.adjust(levene_df$P_value, method = "bonferroni")
      
      incProgress(0.10, detail = "Filtering significant genes...")
      
      # Levene Genes와 Levene Unique 필터링
      levene_fwer <- levene_df %>% 
        filter(`Adjusted P-value` <= input$fwer_value) %>%
        arrange(P_value) %>%
        mutate(`P-value rank` = rank(P_value)) %>%
        select(Gene, `P-value rank`)
      
      levene_unique <- levene_df %>% 
        filter(`Adjusted P-value` <= input$fwer_value & !Gene %in% anova_df$Gene[anova_df$`Adjusted P-value` <= input$fwer_value]) %>%
        arrange(P_value) %>%
        mutate(`P-value rank` = rank(P_value)) %>%
        select(Gene, `P-value rank`)
      
      final_results <- list(
        anova_fwer = anova_df %>% filter(`Adjusted P-value` <= input$fwer_value),
        levene_fwer = levene_fwer,
        anova_unique = anova_df %>% filter(`Adjusted P-value` <= input$fwer_value & !Gene %in% levene_df$Gene[levene_df$`Adjusted P-value` <= input$fwer_value]),
        levene_unique = levene_unique,
        anova_full = anova_df,
        levene_full = levene_df
      )
      
      incProgress(0.05, detail = "Updating results...")
      
      fwer_results(final_results)
      
      updateTabsetPanel(session, "test_result_tabs", selected = "DVGs")
      
      incProgress(0.05, detail = "Complete!")
    })
  })
  
  # Levene Genes 카드 UI
  output$levene_fwer_cards_ui <- renderUI({
    req(fwer_results())
    levene_fwer_data <- fwer_results()$levene_fwer
    if (!all(c("Gene", "P-value rank") %in% colnames(levene_fwer_data))) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; text-align: center; padding: 20px;",
        "No DVGs found."
      ))
    }
    
    # 검색 입력란
    search_input <- tags$input(
      type = "text",
      class = "search-bar",
      placeholder = "Search Genes...",
      oninput = "filterCards(this.value, 'levene-fwer')"
    )
    
    # 카드 리스트
    cards <- lapply(1:nrow(levene_fwer_data), function(i) {
      gene <- levene_fwer_data$Gene[i]
      rank <- levene_fwer_data$`P-value rank`[i]
      div(
        class = "gene-card",
        `data-gene` = tolower(gene),
        onclick = sprintf("Shiny.setInputValue('clicked_gene', '%s')", gene),
        div(class = "gene-card-header", gene),
        div(
          class = "gene-card-content",
          tags$span(sprintf("DVG Rank: %s", rank))
        )
      )
    })
    
    tagList(
      search_input,
      div(class = "cards-container", id = "levene-fwer", cards)
    )
  })
  
  # Levene Unique 카드 UI
  output$levene_unique_cards_ui <- renderUI({
    req(fwer_results())
    levene_unique_data <- fwer_results()$levene_unique
    if (!all(c("Gene", "P-value rank") %in% colnames(levene_unique_data))) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; text-align: center; padding: 20px;",
        "No DVG-only genes found."
      ))
    }
    
    search_input <- tags$input(
      type = "text",
      class = "search-bar",
      placeholder = "Search Genes...",
      oninput = "filterCards(this.value, 'levene-unique')"
    )
    
    cards <- lapply(1:nrow(levene_unique_data), function(i) {
      gene <- levene_unique_data$Gene[i]
      rank <- levene_unique_data$`P-value rank`[i]
      div(
        class = "gene-card",
        `data-gene` = tolower(gene),
        onclick = sprintf("Shiny.setInputValue('clicked_gene', '%s')", gene),
        div(class = "gene-card-header", gene),
        div(
          class = "gene-card-content",
          tags$span(sprintf("DVG-only Rank: %s", rank))
        )
      )
    })
    
    tagList(
      search_input,
      div(class = "cards-container", id = "levene-unique", cards)
    )
  })
  
  
  # 유전자 이미지 및 표 다운로드 핸들러
  output$download_gene_image <- downloadHandler(
    filename = function() {
      paste("gene_analysis_", selected_gene(), "_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      if (is.null(selected_gene()) || is.null(selected_gene_stats()) || 
          is.null(expr_data()) || is.null(coords_data())) {
        showNotification("No gene analysis data available to download.", type = "error")
        return(NULL)
      }
      
      gene <- selected_gene()
      expr <- expr_data()
      coords <- coords_data()
      stats <- selected_gene_stats()
      
      fwer_res <- fwer_results()
      p_value <- NA
      adj_p_value <- NA
      
      if (!is.null(fwer_res) && !is.null(fwer_res$levene_full)) {
        levene_full <- fwer_res$levene_full
        if (gene %in% levene_full$Gene) {
          p_value <- levene_full$P_value[levene_full$Gene == gene]
          adj_p_value <- levene_full$`Adjusted P-value`[levene_full$Gene == gene]
        }
      }
      
      # 과학적 표기법 형식화 함수 (상첨자 표현식용 리스트 반환)
      format_scientific <- function(value) {
        if (length(value) == 0 || is.na(value)) return(NULL)
        sci <- format(value, digits = 5, scientific = TRUE)
        parts <- unlist(strsplit(sci, "e"))
        mantissa <- parts[1]
        exponent <- as.numeric(parts[2])
        list(mantissa = mantissa, exponent = exponent)
      }
      
      p_value_formatted <- format_scientific(p_value)
      adj_p_value_formatted <- format_scientific(adj_p_value)
      
      # 클러스터 데이터
      clusters <- clusters_full()
      cluster_levels <- levels(clusters)
      non_zero_levels <- cluster_levels[cluster_levels != "0"]
      n_colors <- length(non_zero_levels)
      colors <- hcl.colors(n_colors, palette = "Set2")
      color_mapping <- setNames(c("#bdc3c7", colors), c("0", non_zero_levels))
      
      # 플롯 1: 클러스터 분포
      p1 <- ggplot(coords, aes(x = X, y = Y)) +
        geom_point(data = coords[clusters == "0", ], 
                   color = "#bdc3c7", shape = 1, size = 5 * input$test_point_size, alpha = input$test_point_alpha) +
        geom_point(data = coords[clusters != "0", ], 
                   aes(color = factor(clusters[clusters != "0"])), 
                   shape = 16, size = 5 * input$test_point_size, alpha = input$test_point_alpha) +
        scale_color_manual(values = color_mapping) +
        labs(title = "Tissue Domains", x = "X Coordinate", y = "Y Coordinate") +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid = element_blank()
        ) +
        coord_fixed()
      
      # 플롯 2: 유전자 발현
      p2 <- ggplot(coords, aes(x = X, y = Y)) +
        geom_point(data = coords[clusters == "0", ], 
                   color = "#bdc3c7", shape = 1, size = 5 * input$test_point_size, alpha = input$test_point_alpha) +
        geom_point(data = coords[clusters != "0", ], 
                   aes(color = expr[gene, clusters != "0"]), 
                   shape = 16, size = 5 * input$test_point_size, alpha = input$test_point_alpha) +
        scale_color_viridis_c(option = "viridis", name = "Expression") +
        labs(title = paste("Expression of", gene), x = "X Coordinate", y = "Y Coordinate") +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid = element_blank()
        ) +
        coord_fixed()
      
      # 플롯 결합
      plot_row <- plot_grid(p1, p2, ncol = 2, align = "hv", axis = "tblr")
      plot_height <- 300  # 플롯 높이 (포인트 단위)
      plot_width <- 600   # 플롯 너비 (포인트 단위)
      
      # 클러스터 범례
      items <- lapply(1:n_colors, function(i) {
        list(color = colors[i], label = non_zero_levels[i])
      })
      rows <- split(items, ceiling(seq_along(items) / 5))
      legend_height_cluster <- length(rows) * 30  # 각 행당 30pt
      legend_width_cluster <- 300  # 클러스터 범례 너비
      
      png(tempfile())
      pushViewport(viewport(layout = grid.layout(nrow = length(rows), ncol = 5, widths = unit(0.2, "npc"))))
      for (i in seq_along(items)) {
        row <- ceiling(i / 5)
        col <- (i - 1) %% 5 + 1
        pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
        grid.circle(x = 0.2, y = 0.5, r = 0.1, gp = gpar(fill = items[[i]]$color, col = NA))
        grid.text(label = items[[i]]$label, x = 0.4, y = 0.5, just = "left", gp = gpar(fontsize = 12))
        popViewport()
      }
      legend_grob_cluster <- grid.grab()
      dev.off()
      
      # 발현 범례 (가로 방향, "Expression Level" 왼쪽 정렬 유지)
      expr_range <- range(expr[gene, clusters != "0"], na.rm = TRUE)
      min_val <- expr_range[1]
      max_val <- expr_range[2]
      legend_height_expr <- 50  # 발현 범례 높이
      legend_width_expr <- 150  # 발현 범례 너비
      
      png(tempfile())
      gradient_colors <- viridis::viridis(256, option = "viridis")
      gradient_matrix <- matrix(gradient_colors, nrow = 1, ncol = 256)
      pushViewport(viewport(width = unit(150, "points"), height = unit(12, "points")))
      grid.raster(gradient_matrix, width = 1, height = 1, interpolate = TRUE)
      gradient_grob <- grid.grab()
      dev.off()
      
      png(tempfile())
      pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, 
                                                 heights = unit.c(unit(20, "points"), unit(30, "points")))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.text("Expression Level", x = 0, y = 0.5, just = "left", gp = gpar(fontsize = 12, col = "#2c3e50"))
      popViewport()
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, 
                                                 heights = unit.c(unit(12, "points"), unit(18, "points")))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.draw(gradient_grob)
      popViewport()
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      grid.text(sprintf("%.2f", min_val), x = 0, y = 0.5, just = "left", gp = gpar(fontsize = 12, col = "#2c3e50"))
      grid.text(sprintf("%.2f", max_val), x = 1, y = 0.5, just = "right", gp = gpar(fontsize = 12, col = "#2c3e50"))
      popViewport()
      popViewport()
      expr_legend_grob <- grid.grab()
      dev.off()
      
      # 텍스트 영역 크기 계산
      text_height <- 100  # 텍스트 영역 높이 (제목 및 p-value 등)
      text_width <- 600   # 텍스트 영역 너비
      
      # 전체 크기 계산 (표 제외)
      total_width <- max(plot_width, text_width, legend_width_cluster + legend_width_expr + 20)
      legend_height <- max(legend_height_cluster, legend_height_expr) + 20  # 범례 높이 + 여유
      total_height <- text_height + plot_height + legend_height + 30  # 표 제외, 여유 공간 조정
      
      # 사용자 지정 파일 이름 생성
      date_str <- format(Sys.Date(), "%Y%m%d")
      image_file <- paste0("gene_image_", gene, "_", date_str, ".png")
      excel_file <- paste0("gene_stats_", gene, "_", date_str, ".xlsx")
      
      # 임시 파일 경로 생성
      image_path <- file.path(tempdir(), image_file)
      excel_path <- file.path(tempdir(), excel_file)
      
      # 이미지 생성 (표 제외)
      CairoPNG(image_path, width = total_width, height = total_height, units = "pt", dpi = 300)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, 
                                                 heights = unit.c(
                                                   unit(text_height, "points"),
                                                   unit(plot_height, "points"),
                                                   unit(legend_height, "points")
                                                 ))))
      
      # 1. 텍스트 영역 (상첨자 표현식 사용)
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.text(paste("Gene for Analysis:", gene), x = 0.1, y = 0.9, just = "left", 
                gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50"))
      if (!is.null(p_value_formatted)) {
        grid.text(label = bquote("Brown-Forsythe p-value: " ~ .(p_value_formatted$mantissa) ~ " x " ~ 10^.(p_value_formatted$exponent)), 
                  x = 0.1, y = 0.6, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      } else {
        grid.text("p-value: NA", x = 0.1, y = 0.6, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      }
      if (!is.null(adj_p_value_formatted)) {
        grid.text(label = bquote("Adjusted p-value: " ~ .(adj_p_value_formatted$mantissa) ~ " x " ~ 10^.(adj_p_value_formatted$exponent)), 
                  x = 0.1, y = 0.3, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      } else {
        grid.text("Adjusted p-value: NA", x = 0.1, y = 0.3, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      }
      popViewport()
      
      # 2. 플롯 영역
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      grid.draw(ggplotGrob(plot_row))
      popViewport()
      
      # 3. 범례 영역 (좌우 배치)
      pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
      pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, 
                                                 widths = unit.c(unit(legend_width_cluster, "points"), unit(legend_width_expr, "points")))))
      pushViewport(viewport(layout.pos.col = 1))
      grid.draw(legend_grob_cluster)
      popViewport()
      pushViewport(viewport(layout.pos.col = 2))
      grid.draw(expr_legend_grob)
      popViewport()
      popViewport()
      
      dev.off()
      
      # 엑셀 파일 생성 및 스타일링
      stats_table <- stats %>%
        select(Cluster, Mean, Variance) %>%
        mutate(
          Mean = format(round(Mean, 5), nsmall = 5),
          Variance = format(round(Variance, 5), nsmall = 5)
        ) %>%
        rename(Domain = Cluster)
      wb <- createWorkbook()
      addWorksheet(wb, "Stats")
      writeData(wb, "Stats", stats_table, startRow = 1, startCol = 1)
      
      # 스타일 정의
      header_style <- createStyle(
        fontColour = "#ffffff",
        fgFill = "#2980b9",
        fontSize = 11,
        fontName = "Arial",
        textDecoration = "bold",
        halign = "center",
        valign = "center"
      )
      even_row_style <- createStyle(
        fgFill = "#f8f9fa",
        halign = "center",
        valign = "center"
      )
      odd_row_style <- createStyle(
        fgFill = "#ffffff",
        halign = "center",
        valign = "center"
      )
      
      # 스타일 적용
      addStyle(wb, "Stats", style = header_style, rows = 1, cols = 1:3, gridExpand = TRUE)
      for (i in 2:(nrow(stats_table) + 1)) {
        style <- if (i %% 2 == 0) even_row_style else odd_row_style
        addStyle(wb, "Stats", style = style, rows = i, cols = 1:3, gridExpand = TRUE)
      }
      
      # 열 너비 조정
      setColWidths(wb, "Stats", cols = 1:3, widths = "auto")
      
      # 엑셀 파일 저장
      saveWorkbook(wb, excel_path, overwrite = TRUE)
      
      # ZIP 파일 생성
      zip::zipr(file, files = c(image_path, excel_path), include_directories = FALSE)
      
      # 임시 파일 삭제
      unlink(image_path)
      unlink(excel_path)
    }
  )
  
  # "Run Tests" 버튼
  output$test_button_ui <- renderUI({
    is_data_matched <- data_matched()
    expr <- expr_data()
    coords <- coords_data()
    files_present <- !is.null(expr) && !is.null(coords)
    
    btn <- actionButton(
      "test_button",
      tags$span(icon("play"), " Test Start!"),
      class = "btn-primary btn-sm",
      width = "100%"
    )
    
    if (!files_present || !is_data_matched) {
      btn <- tagAppendAttributes(btn, disabled = "disabled")
    }
    
    btn
  })
  
  
  # 유전자 클릭 시 독립적인 분석 수행
  observeEvent(input$clicked_gene, {
    gene <- input$clicked_gene
    selected_gene(gene)
    
    withProgress(message = "Analyzing Selected Gene...", value = 0, {
      req(expr_data(), spots_data(), gene)
      expr <- expr_data()
      spots <- spots_data()
      
      # 데이터 크기 일치 확인
      if (ncol(expr) != nrow(spots)) {
        showNotification(
          paste("Expression matrix columns (", ncol(expr), ") do not match spots (", nrow(spots), ")."),
          type = "error",
          duration = 10
        )
        selected_gene_stats(NULL)
        return()
      }
      
      # 클러스터 소스 선택
      if (input$cluster_source == "brush") {
        clusters <- factor(brush_clusters())
      } else if (input$cluster_source == "file") {
        req(uploaded_clusters())
        clusters <- factor(uploaded_clusters())
        if (length(clusters) != ncol(expr)) {
          stop("Uploaded domain label vector length does not match expression matrix columns.")
        }
      }
      
      # 클러스터 0 제외
      selected_spots <- clusters != 0
      if (length(selected_spots) != ncol(expr)) {
        showNotification("Domain label length does not match expression matrix columns.", type = "error")
        selected_gene_stats(NULL)
        return()
      }
      expr_selected <- expr[, selected_spots, drop = FALSE]
      clusters <- factor(clusters[selected_spots])
      
      if (!gene %in% rownames(expr_selected)) {
        selected_gene_stats(NULL)
        showNotification(paste("Gene", gene, "not found in expression data."), type = "warning")
        return()
      }
      
      incProgress(0.5, detail = "Calculating statistics...")
      gene_data <- expr_selected[gene, ]
      mean_vals <- tapply(gene_data, clusters, mean, na.rm = TRUE)
      var_vals <- tapply(gene_data, clusters, var, na.rm = TRUE)
      df <- data.frame(
        Cluster = names(mean_vals),
        Mean = as.numeric(mean_vals),
        Variance = as.numeric(var_vals),
        Gene = gene
      )
      
      selected_gene_stats(df)
      incProgress(1.0, detail = "Complete")
    })
  })
  
  # 선택된 유전자 분석 UI
  output$test_result_gene_analysis <- renderUI({
    gene <- selected_gene()
    if (is.null(gene)) {
      return(div(
        class = "no-data",
        style = "text-align: center; color: #7f8c8d; padding: 20px;",
        "No gene selected. Click a gene from the table on the left to view details."
      ))
    }
    
    gene_stats <- selected_gene_stats()
    if (is.null(gene_stats)) {
      return(div(
        class = "no-data",
        style = "text-align: center; color: #e74c3c; padding: 20px;",
        paste("No data available for gene:", gene)
      ))
    }
    
    fwer_res <- fwer_results()
    p_value <- NA
    adj_p_value <- NA
    if (!is.null(fwer_res$levene_full)) {
      levene_full <- fwer_res$levene_full
      if (gene %in% levene_full$Gene) {
        p_value <- levene_full$P_value[levene_full$Gene == gene]
        adj_p_value <- levene_full$`Adjusted P-value`[levene_full$Gene == gene]
      }
    }
    
    p_value_text <- format_scientific_notation(p_value)
    adj_p_value_text <- format_scientific_notation(adj_p_value)
    
    tagList(
      div(
        class = "gene-analysis-container",
        h3(
          tags$span("Gene for Analysis: "),
          tags$a(
            href = sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene),
            target = "_blank",
            class = "clickable-gene",
            gene
          ),
          style = "font-size: 18px; color: #2c3e50; margin-bottom: 5px;"
        ),
        p(paste("Brown-Forsythe p-value:", p_value_text), style = "font-size: 14px; color: #2c3e50; margin-bottom: 2px;"),
        p(paste("Adjusted p-value:", adj_p_value_text), style = "font-size: 14px; color: #2c3e50; margin-bottom: 15px;"),
        fluidRow(
          column(6, 
                 plotOutput("test_result_spatial_cluster_plot", height = "250px"),
                 style = "padding: 10px; margin-bottom: 0;"
          ),
          column(6, 
                 plotOutput("test_result_expr_plot", height = "250px"),
                 style = "padding: 10px; margin-bottom: 0;"
          )
        ),
        fluidRow(
          column(6, 
                 uiOutput("test_result_cluster_legend"),
                 style = "text-align: left; padding: 10px;"
          ),
          column(6, 
                 uiOutput("test_result_expr_legend"),
                 style = "text-align: right; padding: 10px;"
          )
        ),
        div(
          style = "margin: 0; padding: 0; max-height: 300px; overflow-y: auto; margin-top: 20px;",
          DTOutput("test_result_stats_table", width = "100%")
        ),
        uiOutput("download_gene_image_ui")
      )
    )
  })
  
  # 반응형 clusters_full 정의
  clusters_full <- reactive({
    if (input$cluster_source == "brush") {
      factor(brush_clusters(), levels = as.character(sort(unique(brush_clusters()))))
    } else if (input$cluster_source == "file") {
      req(uploaded_clusters())
      factor(uploaded_clusters(), levels = as.character(sort(unique(uploaded_clusters()))))
    }
  })
  
  # 선택된 유전자 플롯과 테이블 출력
  output$test_result_spatial_cluster_plot <- renderPlot({
    req(expr_data(), coords_data(), selected_gene(), selected_gene_stats())
    clusters <- clusters_full()
    create_cluster_plot(coords_data(), clusters, input$test_point_size, input$test_point_alpha)
  }, height = 250)
  
  output$test_result_expr_plot <- renderPlot({
    req(expr_data(), coords_data(), selected_gene(), selected_gene_stats())
    gene <- selected_gene()
    expr <- expr_data()
    clusters <- clusters_full()
    create_expression_plot(coords_data(), clusters, as.numeric(expr[gene, ]), gene, 
                           input$test_point_size, input$test_point_alpha)
  }, height = 250)
  
  output$test_result_cluster_legend <- renderUI({
    req(expr_data(), coords_data(), selected_gene(), selected_gene_stats())
    create_cluster_legend(clusters_full())
  })
  
  output$test_result_expr_legend <- renderUI({
    req(expr_data(), coords_data(), selected_gene(), selected_gene_stats())
    gene <- selected_gene()
    expr <- expr_data()
    clusters <- clusters_full()
    df <- data.frame(Expression = as.numeric(expr[gene, ]), Cluster = clusters)
    expr_range <- range(df$Expression[clusters != "0"], na.rm = TRUE)
    create_expression_legend(expr_range)
  })
  
  output$test_result_stats_table <- renderDT({
    req(selected_gene_stats())
    df <- selected_gene_stats() %>% 
      select(Cluster, Mean, Variance) %>%
      mutate(
        Mean = format(round(Mean, 5), nsmall = 5),
        Variance = format(round(Variance, 5), nsmall = 5)
      ) %>%
      rename(Domain = Cluster)
    datatable(
      df,
      options = list(
        dom = 't',              # 테이블만 표시 (검색, 필터 등 제거)
        paging = FALSE,         # 페이지네이션 비활성화
        scrollY = "200px",      # 세로 스크롤 활성화 (높이 조정 가능)
        scrollCollapse = TRUE,  # 스크롤 영역 축소 방지
        autoWidth = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames = FALSE,
      class = "compact"
    )
  })
  
  # Levene Genes 탭 다운로드 UI
  output$download_levene_fwer_ui <- renderUI({
    tagList(
      div(
        style = "clear: both; margin-top: 20px; margin-bottom: 20px;", # 상하단 마진 추가 및 float 해제
        uiOutput("download_results_ui")
      )
    )
  })
  
  # Levene Unique 탭 다운로드 UI
  output$download_levene_unique_ui <- renderUI({
    tagList(
      div(
        style = "clear: both; margin-top: 20px; margin-bottom: 20px;",
        uiOutput("download_results_ui")
      )
    )
  })
  
  # 유전자 이미지 다운로드 버튼 UI
  output$download_gene_image_ui <- renderUI({
    if (is.null(selected_gene()) || is.null(selected_gene_stats()) || 
        is.null(expr_data()) || is.null(coords_data())) {
      div(
        style = "margin-top: 20px; color: #7f8c8d; font-style: italic;",
        "No gene analysis data available."
      )
    } else {
      downloadButton(
        "download_gene_image",
        "Download Gene Analysis File (Image + Table)",
        class = "btn btn-primary btn-sm",
        style = "margin-top: 20px; width: 100%;"
      )
    }
  })
  
  # 통합 다운로드 버튼 UI
  output$download_results_ui <- renderUI({
    if (is.null(fwer_results())) {
      div(
        style = "color: #7f8c8d; font-style: italic; text-align: center;",
        "No test results available."
      )
    } else {
      label <- if (input$test_result_tabs == "DVGs") {
        "Download DVGs Results"
      } else {
        "Download DVG-only Results"
      }
      div(
        style = "text-align: center;", # 버튼 중앙 정렬
        downloadButton(
          "download_results",
          label,
          class = "btn btn-primary btn-sm",
          style = "width: 100%; max-width: 300px;" # 버튼 너비 제한
        )
      )
    }
  })
  
  # 통합 다운로드 핸들러
  output$download_results <- downloadHandler(
    filename = function() {
      if (input$test_result_tabs == "DVGs") {
        paste("DVGs_results_", Sys.Date(), ".csv", sep = "")
      } else {
        paste("DVG_only_results_", Sys.Date(), ".csv", sep = "")
      }
    },
    content = function(file) {
      req(fwer_results())
      levene_full <- fwer_results()$levene_full
      
      if (input$test_result_tabs == "DVGs") {
        levene_fwer <- fwer_results()$levene_fwer
        data <- levene_full %>%
          filter(Gene %in% levene_fwer$Gene) %>%
          mutate(DVG_only = Gene %in% fwer_results()$levene_unique$Gene) %>%
          select(Gene, P_value, `Adjusted P-value`, DVG_only) %>%
          arrange(P_value)
      } else {
        levene_unique <- fwer_results()$levene_unique
        data <- levene_full %>%
          filter(Gene %in% levene_unique$Gene) %>%
          mutate(DVG_only = TRUE) %>%
          select(Gene, P_value, `Adjusted P-value`, DVG_only) %>%
          arrange(P_value)
      }
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  ##############################################################################
  ### STEP 4: DVGs based Tissue Analysis
  
  # Tissue Analysis용 반응형 변수 추가
  tissue_analysis_results <- reactiveVal(NULL)
  tissue_cluster_data <- reactiveVal(NULL)
  
  # Tissue cluster status section (cards only, no title)
  output$tissue_cluster_status_section_cards_only <- renderUI({
    results <- tissue_analysis_results()
    if (is.null(results)) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; padding: 10px; text-align: center;",
        "No spatial domain variability results available."
      ))
    }
    
    cluster_ranking <- results$cluster_ranking %>%
      arrange(Overall_Rank)
    
    # Get cluster colors using viridis
    clusters <- results$clusters_data
    non_zero_clusters <- unique(clusters[clusters != "0"])
    n_colors <- length(non_zero_clusters)
    
    if (n_colors > 0) {
      # Create cluster color mapping based on rank mean
      cluster_ranks <- cluster_ranking$Rank_Mean
      names(cluster_ranks) <- cluster_ranking$Cluster
      
      # Normalize ranks for color mapping
      min_rank <- min(cluster_ranks, na.rm = TRUE)
      max_rank <- max(cluster_ranks, na.rm = TRUE)
      
      if (max_rank > min_rank) {
        normalized_ranks <- (max_rank - cluster_ranks) / (max_rank - min_rank)
      } else {
        normalized_ranks <- rep(0.5, length(cluster_ranks))
        names(normalized_ranks) <- names(cluster_ranks)
      }
      
      cluster_colors <- viridis::viridis(100, option = "viridis")[round(normalized_ranks * 99) + 1]
      names(cluster_colors) <- names(cluster_ranks)
    } else {
      cluster_colors <- character(0)
    }
    
    # Create gradient bar for reference
    gradient_bar <- div(
      style = "margin-bottom: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
      h6("Variability", style = "margin-bottom: 8px; color: #2980b9; font-size: 12px;"),
      div(
        style = "position: relative; width: 100%; height: 15px; background: linear-gradient(to right, #440154, #31688E, #35B779, #FDE725); border-radius: 3px; margin-bottom: 5px;"
      ),
      div(
        style = "display: flex; justify-content: space-between; font-size: 10px; color: #666;",
        tags$span("Low"),
        tags$span("High")
      )
    )
    
    cards <- lapply(1:nrow(cluster_ranking), function(i) {
      cluster <- cluster_ranking$Cluster[i]
      rank <- cluster_ranking$Overall_Rank[i]
      rank_mean <- round(cluster_ranking$Rank_Mean[i], 2)
      iqr_mean <- round(cluster_ranking$IQR_Mean[i], 2)
      
      cluster_color <- if (cluster %in% names(cluster_colors)) {
        cluster_colors[cluster]
      } else {
        "#666666"
      }
      
      div(
        class = "gene-card tissue-cluster-card",
        id = paste0("tissue-cluster-card-", cluster),
        onclick = sprintf("
        var timestamp = new Date().getTime(); 
        Shiny.setInputValue('clicked_tissue_cluster_card', '%s_' + timestamp, {priority: 'event'});
      ", cluster),
        style = paste0("border-left: 4px solid ", cluster_color, "; max-width: 220px; margin-bottom: 8px;"),
        div(
          class = "gene-card-header", 
          style = paste0("color: ", cluster_color, "; display: flex; align-items: center; font-size: 12px;"),
          tags$span(
            style = paste0("display: inline-block; width: 10px; height: 10px; border-radius: 50%; background-color: ", cluster_color, "; margin-right: 6px;")
          ),
          paste("Domain", cluster)
        ),
        div(
          class = "gene-card-content",
          style = "font-size: 10px;",
          tags$span(paste("Rank:", rank)),
          tags$br(),
          tags$span(sprintf("Mean Rank: %.2f", rank_mean))
        )
      )
    })
    
    tagList(
      gradient_bar,
      div(
        class = "cards-container", 
        id = "tissue-cluster-status", 
        style = "max-height: 400px; max-width: 280px; overflow-y: auto;",
        cards
      )
    )
  })
  
  # Tissue cluster card click handler
  observeEvent(input$clicked_tissue_cluster_card, {
    req(tissue_analysis_results(), input$clicked_tissue_cluster_card)
    
    # Extract cluster ID from the clicked card
    cluster_id_raw <- input$clicked_tissue_cluster_card
    cluster_id <- tryCatch({
      cluster_part <- strsplit(as.character(cluster_id_raw), "_")[[1]][1]
      as.numeric(cluster_part)  # 숫자로 변환
    }, error = function(e) {
      showNotification("Invalid domain ID", type = "error")
      return(NULL)
    })
    
    if (is.null(cluster_id) || is.na(cluster_id)) {
      return()
    }
    
    session$sendCustomMessage("highlightTissueCluster", list(
      highlightedCluster = cluster_id  # 숫자로 전달
    ))
  })
  
  # Observer to initialize tissue canvas when analysis is complete
  observe({
    results <- tissue_analysis_results()
    
    if (!is.null(results)) {
      # 스케일된 좌표 데이터 사용
      spots <- results$spots_data
      if (is.null(spots)) {
        spots <- spots_data()  # fallback
      }
      
      if (is.null(spots)) {
        return()
      }
      
      clusters <- results$clusters_data
      cluster_ranking <- results$cluster_ranking
      
      # Prepare spots data for canvas
      rank_means <- cluster_ranking$Rank_Mean
      names(rank_means) <- cluster_ranking$Cluster
      
      # Non-zero clusters and normalization
      non_zero_clusters <- unique(clusters[clusters != "0"])
      if (length(non_zero_clusters) == 0) {
        return()
      }
      
      min_rank <- min(rank_means, na.rm = TRUE)
      max_rank <- max(rank_means, na.rm = TRUE)
      
      if (max_rank > min_rank) {
        normalized_ranks <- (max_rank - rank_means) / (max_rank - min_rank)
      } else {
        normalized_ranks <- rep(0.5, length(rank_means))
        names(normalized_ranks) <- names(rank_means)
      }
      
      # Colors
      cluster_colors <- viridis::viridis(100, option = "viridis")[round(normalized_ranks * 99) + 1]
      
      # named vector를 named list로 변환
      cluster_colors_list <- as.list(cluster_colors)
      names(cluster_colors_list) <- names(normalized_ranks)
      
      # Create spots list for JavaScript - 스케일된 좌표 사용
      spots_list <- lapply(1:nrow(spots), function(i) {
        cluster <- as.character(clusters[i])
        
        if (cluster != "0" && cluster %in% names(cluster_colors_list)) {
          color <- cluster_colors_list[[cluster]]
        } else {
          color <- "#bdc3c7"
        }
        
        list(
          X = spots$X[i],  # 스케일된 X 좌표
          Y = spots$Y[i],  # 스케일된 Y 좌표
          cluster = cluster,
          color = color,
          alpha = 0.8
        )
      })
      
      # Create cluster ranking info for tooltips
      cluster_ranking_info <- list()
      for (i in 1:nrow(cluster_ranking)) {
        cluster <- cluster_ranking$Cluster[i]
        cluster_ranking_info[[cluster]] <- list(
          rank = cluster_ranking$Overall_Rank[i],
          rankMean = round(cluster_ranking$Rank_Mean[i], 3)
        )
      }
      
      # Send data to JavaScript
      session$sendCustomMessage("initTissueCanvas", list(
        spots = spots_list,
        clusterRanking = cluster_ranking_info
      ))
    }
  })
  
  # Tissue Analysis 준비 여부
  output$tissue_analysis_ready <- reactive({
    !is.null(fwer_results()) && nrow(fwer_results()$levene_fwer) > 0
  })
  outputOptions(output, "tissue_analysis_ready", suspendWhenHidden = FALSE)
  
  # Run Tissue Analysis 버튼
  output$run_tissue_analysis_ui <- renderUI({
    is_data_matched <- data_matched()
    files_present <- !is.null(expr_data()) && !is.null(coords_data())
    
    btn <- actionButton(
      "run_tissue_analysis",
      tags$span(icon("play"), " Run Analysis"),
      class = "btn-primary btn-sm",
      width = "100%"
    )
    
    if (!files_present || !is_data_matched) {
      btn <- tagAppendAttributes(btn, disabled = "disabled")
    }
    
    btn
  })
  
  # Tissue Analysis 실행
  observeEvent(input$run_tissue_analysis, {
    if (!check_open_unit_interval(input$posthoc_alpha, "Significance level")) return()
    withProgress(message = "Running Spatial Domain Variability Test...", value = 0, {
      req(data_matched(), fwer_results(), expr_data(), spots_data())
      
      fwer_res <- fwer_results()
      expr <- expr_data()
      spots <- spots_data()
      
      incProgress(0.05, detail = "Preparing DVGs list...")
      
      levene_genes <- fwer_res$levene_fwer %>% 
        arrange(`P-value rank`) %>% 
        pull(Gene)
      
      if (length(levene_genes) == 0) {
        showNotification("No DVGs available.", type = "error")
        return()
      }
      
      incProgress(0.05, detail = "Preparing domain data...")
      
      if (input$cluster_source == "brush") {
        clusters <- factor(brush_clusters())
      } else if (input$cluster_source == "file") {
        req(uploaded_clusters())
        clusters <- factor(uploaded_clusters())
        if (length(clusters) != ncol(expr)) {
          stop("Uploaded domain label vector length does not match expression matrix columns.")
        }
      }
      
      selected_spots <- clusters != 0
      expr_selected <- expr[, selected_spots, drop = FALSE]
      clusters_selected <- factor(clusters[selected_spots])
      
      if (nlevels(clusters_selected) < 2) {
        showNotification("At least 2 domains are required for analysis.", type = "error")
        return()
      }
      
      cluster_levels <- levels(clusters_selected)
      n_genes <- length(levene_genes)
      n_clusters <- length(cluster_levels)
      
      incProgress(0.05, detail = sprintf("Calculating IQRs for %d genes...", n_genes))
      
      # 클러스터별 인덱스 미리 계산
      cluster_indices <- lapply(cluster_levels, function(cl) {
        which(clusters_selected == cl)
      })
      names(cluster_indices) <- cluster_levels
      
      # IQR 행렬 계산 (벡터화)
      iqr_matrix <- matrix(0, nrow = n_genes, ncol = n_clusters)
      rownames(iqr_matrix) <- levene_genes
      colnames(iqr_matrix) <- cluster_levels
      
      for (i in seq_along(levene_genes)) {
        gene <- levene_genes[i]
        if (gene %in% rownames(expr_selected)) {
          gene_data <- expr_selected[gene, ]
          for (j in seq_along(cluster_levels)) {
            cluster_data <- gene_data[cluster_indices[[j]]]
            if (length(cluster_data) > 1) {
              iqr_matrix[i, j] <- IQR(cluster_data, na.rm = TRUE)
            }
          }
        }
      }
      
      incProgress(0.50, detail = "Converting IQRs to ranks...")
      
      rank_matrix <- t(apply(iqr_matrix, 1, function(x) {
        rank(-x, ties.method = "average")
      }))
      
      incProgress(0.10, detail = "Performing Friedman Test...")
      
      friedman_result <- tryCatch({
        friedman.test(rank_matrix)
      }, error = function(e) {
        showNotification(paste("Friedman test failed:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(friedman_result)) return()
      
      posthoc_result <- NULL
      alpha <- input$posthoc_alpha
      if (friedman_result$p.value < alpha) {
        incProgress(0.05, detail = "Performing Nemenyi post-hoc analysis...")
        
        posthoc_result <- tryCatch({
          cluster_ranks <- colMeans(rank_matrix, na.rm = TRUE)
          n_genes_test <- nrow(rank_matrix)
          n_clusters_test <- ncol(rank_matrix)
          q_alpha <- qtukey(1 - alpha, n_clusters_test, Inf)
          critical_diff <- q_alpha * sqrt((n_clusters_test * (n_clusters_test + 1)) / (12 * n_genes_test))
          
          cluster_names <- colnames(rank_matrix)
          comparison_results <- data.frame(
            Cluster1 = character(),
            Cluster2 = character(),
            Mean_Rank_Diff = numeric(),
            Critical_Diff = numeric(),
            Significant = logical(),
            stringsAsFactors = FALSE
          )
          
          for (i in 1:(n_clusters_test-1)) {
            for (j in (i+1):n_clusters_test) {
              cluster1 <- cluster_names[i]
              cluster2 <- cluster_names[j]
              rank_diff <- abs(cluster_ranks[i] - cluster_ranks[j])
              is_significant <- rank_diff > critical_diff
              
              comparison_results <- rbind(comparison_results, data.frame(
                Cluster1 = cluster1,
                Cluster2 = cluster2,
                Mean_Rank_Diff = round(rank_diff, 3),
                Critical_Diff = round(critical_diff, 3),
                Significant = is_significant,
                stringsAsFactors = FALSE
              ))
            }
          }
          
          comparison_results <- comparison_results[order(comparison_results$Mean_Rank_Diff, decreasing = TRUE), ]
          
          list(
            method = "Nemenyi test",
            comparisons = comparison_results,
            cluster_ranks = cluster_ranks,
            critical_difference = critical_diff,
            q_value = q_alpha
          )
        }, error = function(e) {
          showNotification(paste("Nemenyi test failed:", e$message), type = "warning")
          NULL
        })
      } else {
        incProgress(0.05, detail = sprintf("Skipping post-hoc (p >= %s)...", alpha))
      }
      
      incProgress(0.10, detail = "Calculating domain rankings...")
      
      iqr_sums <- colSums(iqr_matrix, na.rm = TRUE)
      iqr_means <- iqr_sums / nrow(iqr_matrix)
      rank_sums <- colSums(rank_matrix, na.rm = TRUE)
      rank_means <- rank_sums / nrow(rank_matrix)
      
      cluster_ranking <- data.frame(
        Cluster = names(iqr_sums),
        IQR_Sum = round(iqr_sums, 3),
        IQR_Mean = round(iqr_means, 3),
        Rank_Sum = rank_sums,
        Rank_Mean = round(rank_means, 3),
        Overall_Rank = rank(rank_means),
        stringsAsFactors = FALSE
      ) %>%
        arrange(Overall_Rank)
      
      incProgress(0.05, detail = "Saving results...")
      
      tissue_results <- list(
        friedman_test = friedman_result,
        posthoc_test = posthoc_result,
        iqr_matrix = iqr_matrix,
        rank_matrix = rank_matrix,
        cluster_ranking = cluster_ranking,
        levene_genes = levene_genes,
        clusters_data = clusters,
        spots_data = spots,
        alpha = alpha
      )
      
      tissue_analysis_results(tissue_results)
      tissue_cluster_data(clusters)
      
      incProgress(0.05, detail = "Complete!")
    })
  })
  
  # Friedman Test 결과 출력
  output$tissue_friedman_result_sidebar <- renderUI({
    results <- tissue_analysis_results()
    if (is.null(results)) return(NULL)
    
    friedman_p <- results$friedman_test$p.value
    chi_squared <- results$friedman_test$statistic
    df <- results$friedman_test$parameter
    
    div(
      style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-top: 10px; border-left: 4px solid #2980b9;",
      h5("Friedman Test Results", style = "color: #2980b9; margin-bottom: 10px; font-size: 14px;"),
      p(sprintf("Chi-squared: %.4f", chi_squared), style = "margin-bottom: 3px; font-size: 12px;"),
      p(sprintf("df: %d", df), style = "margin-bottom: 3px; font-size: 12px;"),
      p(sprintf("p-value: %.2e", friedman_p), 
        style = if (friedman_p < results$alpha) "margin-bottom: 0px; font-weight: bold; color: #d32f2f; font-size: 12px;" else "margin-bottom: 0px; font-weight: bold; font-size: 12px;")
    )
  })
  
  
  # Critical Difference Diagram
  output$tissue_critical_diff_diagram <- renderPlot({
    results <- tissue_analysis_results()
    if (is.null(results) || is.null(results$posthoc_test)) return(NULL)
    
    posthoc <- results$posthoc_test
    if (posthoc$method != "Nemenyi test") return(NULL)
    
    # 클러스터 평균 랭크 정렬 (오름차순 - 가장 좋은 랭크가 왼쪽)
    cluster_ranks <- sort(posthoc$cluster_ranks)
    cluster_names <- names(cluster_ranks)
    n_clusters <- length(cluster_names)
    
    critical_diff <- posthoc$critical_difference
    
    # 수평선의 범위와 중심 계산
    line_start <- min(cluster_ranks) - 0.18
    line_end <- max(cluster_ranks) + 0.18
    line_center <- (line_start + line_end) / 2
    
    # 수평선의 y 위치
    main_line_y <- 0.5
    
    # 라벨을 왼쪽 끝과 오른쪽 끝에 배치
    left_clusters <- cluster_names[cluster_ranks < line_center]
    right_clusters <- cluster_names[cluster_ranks >= line_center]
    
    # 클러스터 개수에 비례한 동적 간격 및 길이 조정
    base_spacing <- max(0.06, min(0.12, 0.05 + n_clusters * 0.009))
    
    # 꺾은선의 시작 높이를 더 가깝게 조정
    start_bend_y <- main_line_y - (0.16 + n_clusters * 0.018)
    
    # 각 클러스터별로 다른 bend_y 높이 설정 (더 작은 간격)
    all_bend_y <- seq(start_bend_y, 
                      start_bend_y - (n_clusters - 1) * base_spacing, 
                      by = -base_spacing)
    
    # 라벨 위치 데이터프레임 생성
    label_positions <- data.frame(
      cluster = character(),
      rank = numeric(),
      label_x = numeric(),
      bend_y = numeric(),
      stringsAsFactors = FALSE
    )
    
    # 왼쪽 클러스터들
    if (length(left_clusters) > 0) {
      left_indices <- which(cluster_names %in% left_clusters)
      left_data <- data.frame(
        cluster = left_clusters,
        rank = cluster_ranks[left_clusters],
        label_x = rep(line_start, length(left_clusters)),
        bend_y = all_bend_y[left_indices],
        stringsAsFactors = FALSE
      )
      label_positions <- rbind(label_positions, left_data)
    }
    
    # 오른쪽 클러스터들
    if (length(right_clusters) > 0) {
      right_indices <- which(cluster_names %in% right_clusters)
      right_data <- data.frame(
        cluster = right_clusters,
        rank = cluster_ranks[right_clusters],
        label_x = rep(line_end, length(right_clusters)),
        bend_y = all_bend_y[right_indices],
        stringsAsFactors = FALSE
      )
      label_positions <- rbind(label_positions, right_data)
    }
    
    # 연결할 클러스터 그룹 찾기
    # mean rank가 낮은 순서대로 정렬 (이미 cluster_ranks가 정렬되어 있음)
    # 먼저 모든 CD 이내의 연결 쌍을 수집
    all_connections <- list()
    
    for (i in 1:n_clusters) {
      current_cluster <- cluster_names[i]
      current_rank <- cluster_ranks[i]
      
      farthest_connected_cluster <- NULL
      
      for (j in (i+1):n_clusters) {
        if (j > n_clusters) break
        
        next_cluster <- cluster_names[j]
        next_rank <- cluster_ranks[j]
        
        if (abs(current_rank - next_rank) <= critical_diff) {
          farthest_connected_cluster <- next_cluster
        } else {
          break
        }
      }
      
      if (!is.null(farthest_connected_cluster)) {
        connected_group <- c(current_cluster, farthest_connected_cluster)
        all_connections[[length(all_connections) + 1]] <- list(
          start = i,
          end = which(cluster_names == farthest_connected_cluster),
          clusters = connected_group
        )
      }
    }
    
    # 포함 관계 필터링: 다른 선에 완전히 포함되는 선 제거
    groups <- list()
    if (length(all_connections) > 0) {
      is_contained <- rep(FALSE, length(all_connections))
      
      for (i in seq_along(all_connections)) {
        conn_i <- all_connections[[i]]
        for (j in seq_along(all_connections)) {
          if (i == j) next
          conn_j <- all_connections[[j]]
          
          # conn_i가 conn_j에 완전히 포함되는지 확인
          if (conn_j$start <= conn_i$start && conn_j$end >= conn_i$end) {
            if (!(conn_j$start == conn_i$start && conn_j$end == conn_i$end)) {
              is_contained[i] <- TRUE
              break
            }
          }
        }
      }
      
      # 포함되지 않는 연결만 groups에 추가
      for (i in seq_along(all_connections)) {
        if (!is_contained[i]) {
          groups[[length(groups) + 1]] <- all_connections[[i]]$clusters
        }
      }
    }
    
    # Y축 범위를 조정 (아래쪽 간격 더 축소)
    min_y_limit <- min(all_bend_y) - 0.18 - length(groups) * 0.035
    max_y_limit <- main_line_y + 0.28  # CD scale 거리는 유지
    
    # 플롯 생성
    p <- ggplot() +
      # 메인 수평선
      geom_segment(aes(x = line_start, xend = line_end,
                       y = main_line_y, yend = main_line_y),
                   color = "#2c3e50", size = 1.2) +
      
      # 각 클러스터에서 수직선 (메인선에서 각자의 bend_y까지)
      geom_segment(data = label_positions,
                   aes(x = rank, xend = rank, y = main_line_y, yend = bend_y),
                   color = "#2c3e50", size = 0.8) +
      
      # 수평선 부분 (각자의 bend_y에서 라벨 x 위치까지)
      geom_segment(data = label_positions,
                   aes(x = rank, xend = label_x, y = bend_y, yend = bend_y),
                   color = "#2c3e50", size = 0.8) +
      
      # 클러스터 라벨 (꺾은선의 끝에 표시)
      geom_text(data = label_positions,
                aes(x = label_x, y = bend_y, label = cluster),
                vjust = 0.5, hjust = ifelse(label_positions$label_x == line_start, 1.1, -0.1), 
                size = 3.5, fontface = "bold") +
      
      # 메인 라인 위의 점들
      geom_point(aes(x = cluster_ranks, y = rep(main_line_y, n_clusters)),
                 color = "#2980b9", size = 3) +
      
      # x축 설정
      scale_x_continuous(limits = c(line_start - 0.18, line_end + 0.18)) +
      
      # 더 컴팩트한 Y축 범위 적용
      scale_y_continuous(limits = c(min_y_limit, max_y_limit)) +
      
      labs(
        title = "Critical Difference Diagram",
        x = "Mean Rank"
      ) +
      
      theme_void() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                  margin = margin(b = 0)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#e74c3c",
                                     margin = margin(b = 2)),
        plot.margin = margin(6, 10, 10, 10)
      )
    
    # Critical Difference 척도 선 추가
    scale_x_start <- line_start
    scale_y <- main_line_y + 0.12
    
    # Critical difference 길이만큼의 수평선
    p <- p + geom_segment(
      aes(x = scale_x_start, xend = scale_x_start + critical_diff, 
          y = scale_y, yend = scale_y),
      color = "#2c3e50", size = 1
    )
    
    # 양쪽 끝의 세로 눈금
    p <- p + geom_segment(
      aes(x = scale_x_start, xend = scale_x_start, 
          y = scale_y - 0.015, yend = scale_y + 0.015),
      color = "#2c3e50", size = 1
    )
    
    p <- p + geom_segment(
      aes(x = scale_x_start + critical_diff, xend = scale_x_start + critical_diff, 
          y = scale_y - 0.015, yend = scale_y + 0.015),
      color = "#2c3e50", size = 1
    )
    
    # 척도 라벨 (적당한 거리로)
    p <- p + geom_text(
      aes(x = scale_x_start + critical_diff/2, y = scale_y + 0.04,
          label = paste("CD =", round(critical_diff, 4))),
      size = 3, hjust = 0.5, color = "#e74c3c"
    )
    
    # Critical difference보다 작은 차이를 가진 그룹들에 대해 수평선 연결
    if (length(groups) > 0) {
      segment_data <- data.frame()
      
      # 그룹 연결선의 간격도 더 작게 조정
      group_spacing <- max(0.03, min(0.05, 0.018 + length(groups) * 0.006))
      
      for (i in seq_along(groups)) {
        group <- groups[[i]]
        if (length(group) > 1) {
          group_indices <- which(cluster_names %in% group)
          group_ranks <- cluster_ranks[group_indices]
          
          if (length(group_ranks) > 0) {
            min_rank <- min(group_ranks)
            max_rank <- max(group_ranks)
            
            connect_y <- main_line_y - 0.05 - (i - 1) * group_spacing
            
            # 데이터프레임에 추가
            segment_data <- rbind(segment_data, data.frame(
              x = min_rank,
              xend = max_rank,
              y = connect_y,
              yend = connect_y,
              group_id = i
            ))
          }
        }
      }
      
      # 모든 빨간선을 한 번에 그리기
      if (nrow(segment_data) > 0) {
        p <- p + geom_segment(
          data = segment_data,
          aes(x = x, xend = xend, y = y, yend = yend),
          color = "#e74c3c", size = 1, alpha = 1
        )
      }
    }
    
    p
  }, height = function() {
    # 동적 높이 계산
    results <- tissue_analysis_results()
    if (is.null(results) || is.null(results$posthoc_test)) return(350)
    
    posthoc <- results$posthoc_test
    if (posthoc$method != "Nemenyi test") return(350)
    
    cluster_ranks <- sort(posthoc$cluster_ranks)
    n_clusters <- length(cluster_ranks)
    critical_diff <- posthoc$critical_difference
    
    # 그룹 개수 계산
    cluster_names <- names(cluster_ranks)
    groups <- list()
    
    for (i in 1:n_clusters) {
      current_cluster <- cluster_names[i]
      current_rank <- cluster_ranks[i]
      
      farthest_connected_cluster <- NULL
      
      for (j in (i+1):n_clusters) {
        if (j > n_clusters) break
        
        next_cluster <- cluster_names[j]
        next_rank <- cluster_ranks[j]
        
        if (abs(current_rank - next_rank) <= critical_diff) {
          farthest_connected_cluster <- next_cluster
        } else {
          break
        }
      }
      
      if (!is.null(farthest_connected_cluster)) {
        connected_group <- c(current_cluster, farthest_connected_cluster)
        groups[[length(groups) + 1]] <- connected_group
      }
    }
    
    # 높이 계산
    base_height <- 280
    cluster_height <- n_clusters * 14
    group_height <- length(groups) * 9
    
    total_height <- base_height + cluster_height + group_height
    
    # 최소/최대 높이
    max(340, min(800, total_height))
  })
  
  
  ##############################################################################
  ### STEP.5 Spatial Expression Pattern Clustering
  
  # Spatial Analysis용 반응형 변수 추가
  levene_graph_data <- reactiveVal(NULL)
  selected_cluster_node <- reactiveVal(NULL)
  selected_levene_gene <- reactiveVal(NULL)
  selected_levene_gene_stats <- reactiveVal(NULL)
  
  # top_n_levene_ui 렌더링
  output$top_n_levene_ui <- renderUI({
    fwer_res <- fwer_results()
    
    # Levene 테스트 결과가 없으면 안내 메시지 표시
    if (is.null(fwer_res) || is.null(fwer_res$levene_fwer) || nrow(fwer_res$levene_fwer) == 0) {
      div(
        style = "color: #7f8c8d; font-style: italic; padding: 10px;",
        "Run STEP 3 (Test for DVGs) first."
      )
    } else {
      # Levene 테스트 결과가 있으면 슬라이더 렌더링
      levene_count <- nrow(fwer_res$levene_fwer)
      step_size <- if (levene_count < 10) {
        1
      } else if (levene_count < 100) {
        5
      } else if (levene_count < 1000) {
        10
      } else {
        50
      }
      max_value <- levene_count
      default_value <- min(100, max_value)
      
      sliderInput(
        "top_n_levene",
        "Number of Top DVGs:",
        min = 0,
        max = max_value,
        value = default_value,
        step = step_size,
        width = "100%",
        ticks = FALSE
      )
    }
  })
  
  # fwer_results 변경 시 슬라이더 상태 관리
  observe({
    fwer_res <- fwer_results()
    
    if (is.null(fwer_res) || is.null(fwer_res$levene_fwer) || nrow(fwer_res$levene_fwer) == 0) {
      # Levene 테스트 결과가 없으면 슬라이더 숨김 (UI에서 처리됨)
      shinyjs::hide("top_n_levene")
    } else {
      # Levene 테스트 결과가 있으면 슬라이더 표시 및 활성화
      shinyjs::show("top_n_levene")
      shinyjs::enable("top_n_levene")
      levene_count <- nrow(fwer_res$levene_fwer)
      step_size <- if (levene_count < 10) {
        1
      } else if (levene_count < 100) {
        5
      } else if (levene_count < 1000) {
        10
      } else {
        50
      }
      max_value <- levene_count
      default_value <- min(100, max_value)
      
      # 슬라이더 속성 업데이트
      updateSliderInput(
        session,
        "top_n_levene",
        min = 0,
        max = max_value,
        value = default_value,
        step = step_size
      )
    }
  })
  
  # Spatial Analysis 준비 여부
  output$spatial_analysis_ready <- reactive({
    !is.null(fwer_results()) && nrow(fwer_results()$levene_fwer) > 0
  })
  outputOptions(output, "spatial_analysis_ready", suspendWhenHidden = FALSE)
  
  # 클러스터링
  observeEvent(input$run_spatial_analysis, {
    withProgress(message = "Running Spatial Analysis...", value = 0, {
      req(data_matched(), fwer_results(), expr_data(), coords_data())
      
      fwer_res <- fwer_results()
      expr <- expr_data()
      
      incProgress(0.02, detail = "Selecting Top DVGs...")
      
      top_n <- input$top_n_levene
      levene_genes <- fwer_res$levene_fwer %>% 
        arrange(`P-value rank`) %>% 
        slice_head(n = top_n) %>% 
        pull(Gene)
      
      if (length(levene_genes) == 0) {
        showNotification("No DVGs available.", type = "error")
        return()
      }
      
      levene_expr <- expr[levene_genes, , drop = FALSE]
      if (nrow(levene_expr) < 2) {
        showNotification("Insufficient genes for clustering.", type = "error")
        return()
      }
      
      incProgress(0.02, detail = "Processing coordinate data...")
      
      coords_raw <- coords_data()
      coords_df <- as.data.frame(coords_raw)
      
      if ("Spot" %in% colnames(coords_df)) {
        rownames(coords_df) <- coords_df$Spot
        coords_df$Spot <- NULL
      } else if (ncol(coords_df) > 2) {
        if (is.character(coords_df[, 1]) || is.factor(coords_df[, 1])) {
          rownames(coords_df) <- coords_df[, 1]
          coords_df <- coords_df[, -1, drop = FALSE]
        }
      }
      
      if (ncol(coords_df) != 2) {
        showNotification("Coordinate data must have exactly 2 columns (X, Y).", type = "error")
        return()
      }
      
      coords_matrix <- matrix(0, nrow = nrow(coords_df), ncol = 2)
      rownames(coords_matrix) <- rownames(coords_df)
      
      for (col_idx in 1:2) {
        col_data <- coords_df[, col_idx]
        if (is.factor(col_data)) col_data <- as.character(col_data)
        coords_matrix[, col_idx] <- as.numeric(col_data)
      }
      
      if (any(is.na(coords_matrix))) {
        valid_rows <- complete.cases(coords_matrix)
        coords_matrix <- coords_matrix[valid_rows, , drop = FALSE]
      }
      
      common_samples <- intersect(colnames(levene_expr), rownames(coords_matrix))
      
      if (length(common_samples) < 10) {
        if (nrow(coords_matrix) == ncol(levene_expr)) {
          sample_names <- paste0("Sample_", 1:nrow(coords_matrix))
          rownames(coords_matrix) <- sample_names
          colnames(levene_expr) <- sample_names
          common_samples <- sample_names
        } else {
          showNotification("Cannot match coordinate and expression data.", type = "error")
          return()
        }
      }
      
      coords_matrix <- coords_matrix[common_samples, , drop = FALSE]
      levene_expr <- levene_expr[, common_samples, drop = FALSE]
      
      n_genes <- nrow(levene_expr)
      
      incProgress(0.03, detail = "Building spatial neighborhoods...")
      
      nn_7 <- nn2(coords_matrix, coords_matrix, k = 7)
      spatial_neighbors_7 <- nn_7$nn.idx
      
      nn_19 <- nn2(coords_matrix, coords_matrix, k = 19)
      spatial_neighbors_19 <- nn_19$nn.idx
      
      nn_37 <- nn2(coords_matrix, coords_matrix, k = 37)
      spatial_neighbors_37 <- nn_37$nn.idx
      
      nn_61 <- nn2(coords_matrix, coords_matrix, k = 61)
      spatial_neighbors_61 <- nn_61$nn.idx
      
      incProgress(0.05, detail = "Computing smoothed features...")
      
      compute_smoothed_features_fast <- function(expr_data, neighbor_idx) {
        n_genes <- nrow(expr_data)
        n_spots <- ncol(expr_data)
        K <- ncol(neighbor_idx)
        
        row_indices <- rep(1:n_spots, times = K)
        col_indices <- as.vector(neighbor_idx)
        
        W_binary <- sparseMatrix(i = row_indices, j = col_indices, x = 1, dims = c(n_spots, n_spots))
        W_mean <- W_binary / K
        W_mean_t <- t(W_mean)
        W_binary_t <- t(W_binary)
        
        expr_mean <- as.matrix(expr_data %*% W_mean_t)
        expr_squared <- expr_data ^ 2
        mean_of_squared <- as.matrix(expr_squared %*% W_mean_t)
        expr_var <- mean_of_squared - (expr_mean ^ 2)
        expr_var[expr_var < 0] <- 0
        expr_var[is.na(expr_var)] <- 0
        
        is_zero <- (expr_data == 0) * 1
        expr_zero_ratio <- as.matrix(is_zero %*% W_mean_t)
        
        thresholds <- rowQuantiles(expr_data, probs = 0.90, na.rm = TRUE)
        is_hotspot <- (expr_data - thresholds) >= 0
        mode(is_hotspot) <- "integer"
        neighbor_hotspot_count <- as.matrix(is_hotspot %*% W_binary_t)
        expr_hotspot <- (neighbor_hotspot_count > 0) * 1
        
        rownames(expr_mean) <- rownames(expr_data)
        colnames(expr_mean) <- colnames(expr_data)
        rownames(expr_var) <- rownames(expr_data)
        colnames(expr_var) <- colnames(expr_data)
        rownames(expr_hotspot) <- rownames(expr_data)
        colnames(expr_hotspot) <- colnames(expr_data)
        rownames(expr_zero_ratio) <- rownames(expr_data)
        colnames(expr_zero_ratio) <- colnames(expr_data)
        
        return(list(mean = expr_mean, var = expr_var, hotspot = expr_hotspot, zero_ratio = expr_zero_ratio))
      }
      
      smoothed_7 <- compute_smoothed_features_fast(levene_expr, spatial_neighbors_7)
      expr_mean_7 <- smoothed_7$mean
      expr_var_7 <- smoothed_7$var
      expr_hotspot_7 <- smoothed_7$hotspot
      expr_zero_ratio_7 <- smoothed_7$zero_ratio
      
      smoothed_19 <- compute_smoothed_features_fast(levene_expr, spatial_neighbors_19)
      expr_mean_19 <- smoothed_19$mean
      expr_var_19 <- smoothed_19$var
      expr_hotspot_19 <- smoothed_19$hotspot
      expr_zero_ratio_19 <- smoothed_19$zero_ratio
      
      smoothed_37 <- compute_smoothed_features_fast(levene_expr, spatial_neighbors_37)
      expr_mean_37 <- smoothed_37$mean
      expr_var_37 <- smoothed_37$var
      expr_hotspot_37 <- smoothed_37$hotspot
      expr_zero_ratio_37 <- smoothed_37$zero_ratio
      
      smoothed_61 <- compute_smoothed_features_fast(levene_expr, spatial_neighbors_61)
      expr_mean_61 <- smoothed_61$mean
      expr_var_61 <- smoothed_61$var
      expr_hotspot_61 <- smoothed_61$hotspot
      expr_zero_ratio_61 <- smoothed_61$zero_ratio
      
      incProgress(0.05, detail = "Preparing original expression features...")
      
      expr_zero_binary <- matrix(0, nrow = nrow(levene_expr), ncol = ncol(levene_expr))
      rownames(expr_zero_binary) <- rownames(levene_expr)
      colnames(expr_zero_binary) <- colnames(levene_expr)
      expr_zero_binary[levene_expr == 0] <- 1
      
      expr_top10 <- matrix(0, nrow = nrow(levene_expr), ncol = ncol(levene_expr))
      rownames(expr_top10) <- rownames(levene_expr)
      colnames(expr_top10) <- colnames(levene_expr)
      
      for (i in 1:nrow(levene_expr)) {
        gene_expr <- levene_expr[i, ]
        threshold_90 <- quantile(gene_expr, probs = 0.90, na.rm = TRUE)
        expr_top10[i, ] <- ifelse(gene_expr >= threshold_90, gene_expr, 0)
      }
      
      K <- ceiling(sqrt(n_genes))
      
      build_knn_filter <- function(expr_data, K) {
        hnsw_idx <- hnsw_build(expr_data, distance = "cosine", M = 32, ef = 200)
        hnsw_res <- hnsw_search(expr_data, hnsw_idx, k = K + 1, ef = max(K * 2, 200))
        knn_idx <- hnsw_res$idx[, -1, drop = FALSE]
        
        n_genes <- nrow(expr_data)
        knn_adj <- matrix(0, nrow = n_genes, ncol = n_genes)
        for (i in 1:n_genes) {
          for (j in knn_idx[i, ]) {
            knn_adj[i, j] <- 1
            knn_adj[j, i] <- 1
          }
        }
        diag(knn_adj) <- 0
        return(knn_adj)
      }
      
      incProgress(0.10, detail = "Building 19 KNN filters...")
      
      knn_adj_mean_7 <- build_knn_filter(expr_mean_7, K)
      knn_adj_var_7 <- build_knn_filter(expr_var_7, K)
      knn_adj_hotspot_7 <- build_knn_filter(expr_hotspot_7, K)
      knn_adj_zero_ratio_7 <- build_knn_filter(expr_zero_ratio_7, K)
      
      knn_adj_mean_19 <- build_knn_filter(expr_mean_19, K)
      knn_adj_var_19 <- build_knn_filter(expr_var_19, K)
      knn_adj_hotspot_19 <- build_knn_filter(expr_hotspot_19, K)
      knn_adj_zero_ratio_19 <- build_knn_filter(expr_zero_ratio_19, K)
      
      knn_adj_mean_37 <- build_knn_filter(expr_mean_37, K)
      knn_adj_var_37 <- build_knn_filter(expr_var_37, K)
      knn_adj_hotspot_37 <- build_knn_filter(expr_hotspot_37, K)
      knn_adj_zero_ratio_37 <- build_knn_filter(expr_zero_ratio_37, K)
      
      knn_adj_mean_61 <- build_knn_filter(expr_mean_61, K)
      knn_adj_var_61 <- build_knn_filter(expr_var_61, K)
      knn_adj_hotspot_61 <- build_knn_filter(expr_hotspot_61, K)
      knn_adj_zero_ratio_61 <- build_knn_filter(expr_zero_ratio_61, K)
      
      knn_adj_zero_binary <- build_knn_filter(expr_zero_binary, K)
      knn_adj_original <- build_knn_filter(levene_expr, K)
      knn_adj_top10 <- build_knn_filter(expr_top10, K)
      
      incProgress(0.05, detail = "Combining KNN filters...")
      
      knn_adj_support <- knn_adj_mean_7 + knn_adj_var_7 + knn_adj_hotspot_7 + knn_adj_zero_ratio_7 +
        knn_adj_mean_19 + knn_adj_var_19 + knn_adj_hotspot_19 + knn_adj_zero_ratio_19 +
        knn_adj_mean_37 + knn_adj_var_37 + knn_adj_hotspot_37 + knn_adj_zero_ratio_37 +
        knn_adj_mean_61 + knn_adj_var_61 + knn_adj_hotspot_61 + knn_adj_zero_ratio_61 +
        knn_adj_zero_binary + knn_adj_original + knn_adj_top10
      
      knn_adj_union <- (knn_adj_support > 0) * 1
      
      if (sum(knn_adj_union) == 0) {
        showNotification("No edges after KNN union.", type = "error")
        return()
      }
      
      incProgress(0.08, detail = "Applying Disparity Filter...")
      
      disparity_filter <- function(weight_matrix, alpha = 0.05) {
        n <- nrow(weight_matrix)
        filtered_adj <- matrix(0, n, n)
        
        for (i in 1:n) {
          neighbors_weight <- weight_matrix[i, ]
          neighbors_weight[i] <- NA
          valid_neighbors <- which(!is.na(neighbors_weight) & neighbors_weight > 0)
          
          if (length(valid_neighbors) == 0) next
          k <- length(valid_neighbors)
          if (k <= 1) {
            filtered_adj[i, valid_neighbors] <- 1
            next
          }
          
          neighbor_weights <- neighbors_weight[valid_neighbors]
          sum_weights <- sum(neighbor_weights, na.rm = TRUE)
          if (sum_weights == 0) next
          
          normalized_weights <- neighbor_weights / sum_weights
          
          for (idx in 1:length(valid_neighbors)) {
            j <- valid_neighbors[idx]
            p_ij <- normalized_weights[idx]
            alpha_ij <- (1 - p_ij)^(k - 1)
            if (alpha_ij < alpha) {
              filtered_adj[i, j] <- 1
            }
          }
        }
        
        filtered_adj <- pmax(filtered_adj, t(filtered_adj))
        return(filtered_adj)
      }
      
      disparity_filtered_adj <- disparity_filter(knn_adj_support, alpha = 0.05)
      n_edges_after <- sum(disparity_filtered_adj) / 2
      
      if (n_edges_after == 0) {
        showNotification("No edges after disparity filter.", type = "error")
        return()
      }
      
      filtered_adj_binary <- disparity_filtered_adj
      
      incProgress(0.02, detail = "Extracting filtered edges...")
      
      edge_indices <- which(filtered_adj_binary == 1 & upper.tri(filtered_adj_binary), arr.ind = TRUE)
      n_filtered_edges <- nrow(edge_indices)
      
      if (n_filtered_edges == 0) {
        showNotification("No edges to compute similarities for.", type = "error")
        return()
      }
      
      incProgress(0.03, detail = "Preparing Jaccard sets (4 measures, 37-neighbor)...")
      
      # Using 37-neighbor smoothed features (excluding Jaccard-Zero)
      hotspot_mean_37_sets <- lapply(1:n_genes, function(i) {
        which(expr_mean_37[i, ] >= quantile(expr_mean_37[i, ], probs = 0.90, na.rm = TRUE))
      })
      names(hotspot_mean_37_sets) <- levene_genes
      
      top10_var_37_sets <- lapply(1:n_genes, function(i) {
        which(expr_var_37[i, ] >= quantile(expr_var_37[i, ], probs = 0.90, na.rm = TRUE))
      })
      names(top10_var_37_sets) <- levene_genes
      
      bottom10_var_37_sets <- lapply(1:n_genes, function(i) {
        which(expr_var_37[i, ] <= quantile(expr_var_37[i, ], probs = 0.10, na.rm = TRUE))
      })
      names(bottom10_var_37_sets) <- levene_genes
      
      incProgress(0.02, detail = sprintf("Computing 4 similarities for %d edges (parallel)...", n_filtered_edges))
      
      # 병렬 처리 설정
      num_cores <- max(1, detectCores() - 1)
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      clusterExport(cl, c("edge_indices", "expr_mean_37", "expr_var_37",
                          "hotspot_mean_37_sets", "top10_var_37_sets", 
                          "bottom10_var_37_sets", "levene_genes"), 
                    envir = environment())
      
      # 4 measures only (excluding Jaccard-Zero), using 37-neighbor
      edge_similarities <- foreach(idx = 1:n_filtered_edges, .combine = 'rbind',
                                   .packages = c('base', 'stats')) %dopar% {
                                     i <- edge_indices[idx, 1]
                                     j <- edge_indices[idx, 2]
                                     
                                     # 1. Pearson-Mean (37-neighbor)
                                     cor_val <- cor(expr_mean_37[i, ], expr_mean_37[j, ], method = "pearson", use = "complete.obs")
                                     pearson_mean_37 <- if (is.na(cor_val)) 0 else cor_val
                                     
                                     # 2. Jaccard-Hotspot (37-neighbor)
                                     setA_hot <- hotspot_mean_37_sets[[i]]
                                     setB_hot <- hotspot_mean_37_sets[[j]]
                                     jaccard_hotspot_mean_37 <- if (length(setA_hot) == 0 || length(setB_hot) == 0) 0 else length(intersect(setA_hot, setB_hot)) / length(union(setA_hot, setB_hot))
                                     
                                     # 3. Jaccard-Var-Top (37-neighbor)
                                     setA_top10 <- top10_var_37_sets[[i]]
                                     setB_top10 <- top10_var_37_sets[[j]]
                                     jaccard_top10_var_37 <- if (length(setA_top10) == 0 || length(setB_top10) == 0) 0 else length(intersect(setA_top10, setB_top10)) / length(union(setA_top10, setB_top10))
                                     
                                     # 4. Jaccard-Var-Bot (37-neighbor)
                                     setA_bottom10 <- bottom10_var_37_sets[[i]]
                                     setB_bottom10 <- bottom10_var_37_sets[[j]]
                                     jaccard_bottom10_var_37 <- if (length(setA_bottom10) == 0 || length(setB_bottom10) == 0) 0 else length(intersect(setA_bottom10, setB_bottom10)) / length(union(setA_bottom10, setB_bottom10))
                                     
                                     c(pearson_mean_37, jaccard_hotspot_mean_37, jaccard_top10_var_37, jaccard_bottom10_var_37)
                                   }
      
      stopCluster(cl)
      
      incProgress(0.15, detail = "Building 4 similarity matrices (37-neighbor)...")
      
      pearson_mean_37_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
      jaccard_hotspot_mean_37_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
      jaccard_top10_var_37_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
      jaccard_bottom10_var_37_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
      
      rownames(pearson_mean_37_matrix) <- levene_genes
      colnames(pearson_mean_37_matrix) <- levene_genes
      rownames(jaccard_hotspot_mean_37_matrix) <- levene_genes
      colnames(jaccard_hotspot_mean_37_matrix) <- levene_genes
      rownames(jaccard_top10_var_37_matrix) <- levene_genes
      colnames(jaccard_top10_var_37_matrix) <- levene_genes
      rownames(jaccard_bottom10_var_37_matrix) <- levene_genes
      colnames(jaccard_bottom10_var_37_matrix) <- levene_genes
      
      diag(pearson_mean_37_matrix) <- 1
      diag(jaccard_hotspot_mean_37_matrix) <- 1
      diag(jaccard_top10_var_37_matrix) <- 1
      diag(jaccard_bottom10_var_37_matrix) <- 1
      
      for (idx in 1:n_filtered_edges) {
        i <- edge_indices[idx, 1]
        j <- edge_indices[idx, 2]
        
        pearson_mean_37_matrix[i, j] <- edge_similarities[idx, 1]
        pearson_mean_37_matrix[j, i] <- edge_similarities[idx, 1]
        jaccard_hotspot_mean_37_matrix[i, j] <- edge_similarities[idx, 2]
        jaccard_hotspot_mean_37_matrix[j, i] <- edge_similarities[idx, 2]
        jaccard_top10_var_37_matrix[i, j] <- edge_similarities[idx, 3]
        jaccard_top10_var_37_matrix[j, i] <- edge_similarities[idx, 3]
        jaccard_bottom10_var_37_matrix[i, j] <- edge_similarities[idx, 4]
        jaccard_bottom10_var_37_matrix[j, i] <- edge_similarities[idx, 4]
      }
      
      incProgress(0.08, detail = "Performing 4 Leiden clusterings (37-neighbor)...")
      
      pearson_mean_37_adj <- pearson_mean_37_matrix
      pearson_mean_37_adj[pearson_mean_37_adj < 0] <- 0
      diag(pearson_mean_37_adj) <- 0
      pearson_mean_37_graph <- graph_from_adjacency_matrix(pearson_mean_37_adj, mode = "undirected", weighted = TRUE)
      pearson_mean_37_membership <- membership(cluster_leiden(pearson_mean_37_graph, objective_function = "modularity", resolution = 1))
      
      jaccard_hotspot_mean_37_adj <- jaccard_hotspot_mean_37_matrix
      jaccard_hotspot_mean_37_adj[jaccard_hotspot_mean_37_adj < 0] <- 0
      diag(jaccard_hotspot_mean_37_adj) <- 0
      jaccard_hotspot_mean_37_graph <- graph_from_adjacency_matrix(jaccard_hotspot_mean_37_adj, mode = "undirected", weighted = TRUE)
      jaccard_hotspot_mean_37_membership <- membership(cluster_leiden(jaccard_hotspot_mean_37_graph, objective_function = "modularity", resolution = 1))
      
      jaccard_top10_var_37_adj <- jaccard_top10_var_37_matrix
      jaccard_top10_var_37_adj[jaccard_top10_var_37_adj < 0] <- 0
      diag(jaccard_top10_var_37_adj) <- 0
      jaccard_top10_var_37_graph <- graph_from_adjacency_matrix(jaccard_top10_var_37_adj, mode = "undirected", weighted = TRUE)
      jaccard_top10_var_37_membership <- membership(cluster_leiden(jaccard_top10_var_37_graph, objective_function = "modularity", resolution = 1))
      
      jaccard_bottom10_var_37_adj <- jaccard_bottom10_var_37_matrix
      jaccard_bottom10_var_37_adj[jaccard_bottom10_var_37_adj < 0] <- 0
      diag(jaccard_bottom10_var_37_adj) <- 0
      jaccard_bottom10_var_37_graph <- graph_from_adjacency_matrix(jaccard_bottom10_var_37_adj, mode = "undirected", weighted = TRUE)
      jaccard_bottom10_var_37_membership <- membership(cluster_leiden(jaccard_bottom10_var_37_graph, objective_function = "modularity", resolution = 1))
      
      incProgress(0.08, detail = "Finding consensus modules (4 measures)...")
      
      consensus_adj <- matrix(0, nrow = n_genes, ncol = n_genes)
      rownames(consensus_adj) <- levene_genes
      colnames(consensus_adj) <- levene_genes
      
      # Consensus: all 4 measures must agree (excluding Jaccard-Zero)
      for (idx in 1:n_filtered_edges) {
        i <- edge_indices[idx, 1]
        j <- edge_indices[idx, 2]
        
        if (pearson_mean_37_membership[i] == pearson_mean_37_membership[j] &
            jaccard_hotspot_mean_37_membership[i] == jaccard_hotspot_mean_37_membership[j] &
            jaccard_top10_var_37_membership[i] == jaccard_top10_var_37_membership[j] &
            jaccard_bottom10_var_37_membership[i] == jaccard_bottom10_var_37_membership[j]) {
          consensus_adj[i, j] <- 1
          consensus_adj[j, i] <- 1
        }
      }
      
      consensus_graph <- graph_from_adjacency_matrix(consensus_adj, mode = "undirected")
      consensus_components <- components(consensus_graph)
      consensus_cluster_labels <- consensus_components$membership
      names(consensus_cluster_labels) <- levene_genes
      
      incProgress(0.05, detail = "Filtering by minimum module size...")
      
      min_size <- input$min_cluster_size
      consensus_cluster_counts <- table(consensus_cluster_labels)
      valid_consensus_clusters <- as.numeric(names(consensus_cluster_counts)[consensus_cluster_counts >= min_size])
      final_valid_nodes <- names(consensus_cluster_labels)[consensus_cluster_labels %in% valid_consensus_clusters]
      
      if (length(final_valid_nodes) == 0) {
        showNotification("No modules meet the minimum size requirement.", type = "error")
        return()
      }
      
      final_graph <- induced_subgraph(consensus_graph, vids = final_valid_nodes)
      final_cluster_labels <- consensus_cluster_labels[final_valid_nodes]
      final_levene_genes <- final_valid_nodes
      
      incProgress(0.05, detail = "Building module graph...")
      
      unique_clusters <- unique(final_cluster_labels)
      n_clusters <- length(unique_clusters)
      cluster_graph <- make_empty_graph(n = n_clusters, directed = FALSE)
      V(cluster_graph)$name <- as.character(unique_clusters)
      
      final_cluster_counts <- table(final_cluster_labels)
      cluster_sizes <- as.numeric(final_cluster_counts[as.character(unique_clusters)])
      
      incProgress(0.02, detail = "Saving results...")
      
      graph_data <- list(
        gene_graph = final_graph,
        clusters = final_cluster_labels,
        levene_genes = final_levene_genes,
        cluster_graph = cluster_graph,
        cluster_sizes = cluster_sizes,
        filtered_adjacency = filtered_adj_binary,
        disparity_filtered_adjacency = disparity_filtered_adj,
        spatial_neighbors_7 = spatial_neighbors_7,
        spatial_neighbors_19 = spatial_neighbors_19,
        spatial_neighbors_37 = spatial_neighbors_37,
        spatial_neighbors_61 = spatial_neighbors_61,
        knn_adjacency_mean_7 = knn_adj_mean_7, knn_adjacency_var_7 = knn_adj_var_7,
        knn_adjacency_hotspot_7 = knn_adj_hotspot_7, knn_adjacency_zero_ratio_7 = knn_adj_zero_ratio_7,
        knn_adjacency_mean_19 = knn_adj_mean_19, knn_adjacency_var_19 = knn_adj_var_19,
        knn_adjacency_hotspot_19 = knn_adj_hotspot_19, knn_adjacency_zero_ratio_19 = knn_adj_zero_ratio_19,
        knn_adjacency_mean_37 = knn_adj_mean_37, knn_adjacency_var_37 = knn_adj_var_37,
        knn_adjacency_hotspot_37 = knn_adj_hotspot_37, knn_adjacency_zero_ratio_37 = knn_adj_zero_ratio_37,
        knn_adjacency_mean_61 = knn_adj_mean_61, knn_adjacency_var_61 = knn_adj_var_61,
        knn_adjacency_hotspot_61 = knn_adj_hotspot_61, knn_adjacency_zero_ratio_61 = knn_adj_zero_ratio_61,
        knn_adjacency_zero_binary = knn_adj_zero_binary, knn_adjacency_original = knn_adj_original,
        knn_adjacency_top10 = knn_adj_top10, knn_adjacency_support = knn_adj_support, knn_adjacency_union = knn_adj_union,
        expr_mean_7 = expr_mean_7, expr_var_7 = expr_var_7, expr_hotspot_7 = expr_hotspot_7, expr_zero_ratio_7 = expr_zero_ratio_7,
        expr_mean_19 = expr_mean_19, expr_var_19 = expr_var_19, expr_hotspot_19 = expr_hotspot_19, expr_zero_ratio_19 = expr_zero_ratio_19,
        expr_mean_37 = expr_mean_37, expr_var_37 = expr_var_37, expr_hotspot_37 = expr_hotspot_37, expr_zero_ratio_37 = expr_zero_ratio_37,
        expr_mean_61 = expr_mean_61, expr_var_61 = expr_var_61, expr_hotspot_61 = expr_hotspot_61, expr_zero_ratio_61 = expr_zero_ratio_61,
        expr_zero_binary = expr_zero_binary, expr_top10 = expr_top10,
        consensus_adjacency = consensus_adj, filtered_edge_list = edge_indices, filtered_edge_similarities = edge_similarities,
        individual_consensus_clusters = list(
          pearson_37_mean = pearson_mean_37_membership,
          jaccard_37_mean_hotspot = jaccard_hotspot_mean_37_membership, 
          jaccard_37_var_top10 = jaccard_top10_var_37_membership,
          jaccard_37_var_bottom10 = jaccard_bottom10_var_37_membership
        ),
        pearson_mean_37_matrix = pearson_mean_37_matrix,
        jaccard_hotspot_mean_37_matrix = jaccard_hotspot_mean_37_matrix, 
        jaccard_top10_var_37_matrix = jaccard_top10_var_37_matrix,
        jaccard_bottom10_var_37_matrix = jaccard_bottom10_var_37_matrix
      )
      
      levene_graph_data(graph_data)
      
      incProgress(0.02, detail = "Complete!")
    })
  })
  
  # 클러스터 그래프 시각화
  output$cluster_network_plot <- renderVisNetwork({
    req(data_matched(), levene_graph_data())
    graph_data <- levene_graph_data()
    g <- graph_data$cluster_graph
    cluster_sizes <- graph_data$cluster_sizes
    
    # 그래프가 비어있는지 확인
    if (vcount(g) == 0) {
      return(visNetwork(data.frame(id = character(0)), data.frame()))
    }
    
    # 최대 노드 크기와 엣지 굵기 설정
    max_node_size <- 50
    max_edge_width <- 5
    
    # 노드 크기 스케일링
    max_cluster_size <- max(cluster_sizes)
    scaled_sizes <- if (max_cluster_size > 0) {
      (cluster_sizes / max_cluster_size) * max_node_size
    } else {
      rep(max_node_size, length(cluster_sizes))
    }
    
    # 노드 데이터프레임
    nodes <- data.frame(
      id = V(g)$name,
      label = paste("Module", V(g)$name),
      size = scaled_sizes,
      color = hcl.colors(length(V(g)), palette = "Viridis"),
      title = paste("Module:", V(g)$name, "<br>Genes:", cluster_sizes)
    )
    
    # 엣지가 있는 경우에만 처리
    if (ecount(g) > 0 && !is.null(E(g)$weight)) {
      # 엣지 데이터프레임 - igraph 함수임을 명시적으로 지정
      edge_df <- igraph::as_data_frame(g, what = "edges")
      
      max_edge_weight <- max(E(g)$weight, 1)
      min_edge_weight <- min(E(g)$weight, 0)
      
      edges <- edge_df %>%
        mutate(
          # 거리가 작을수록 더 굵게 (1 - normalized_distance)
          width = ((1 - E(g)$weight) / (1 - min_edge_weight + 0.001)) * max_edge_width,
          # 거리가 작을수록 더 짧게 (거리 값에 비례하여 length 설정)
          length = 50 + (E(g)$weight * 200),  # 최소 50, 최대 250 길이
          color = list(opacity = 0.4),
          title = sprintf("Edge Weight: %.3f", E(g)$weight)
        )
    } else if (ecount(g) > 0) {
      # weight이 없는 경우
      edge_df <- igraph::as_data_frame(g, what = "edges")
      edges <- edge_df %>%
        mutate(
          width = 2,
          length = 150,
          color = list(opacity = 0.4),
          title = "Edge"
        )
    } else {
      # 엣지가 없는 경우 빈 데이터프레임
      edges <- data.frame(
        from = character(0),
        to = character(0),
        width = numeric(0),
        length = numeric(0)
      )
    }
    
    visNetwork(nodes, edges, height = "500px", width = "100%") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visPhysics(enabled = TRUE) %>%
      visLayout(randomSeed = 123, improvedLayout = TRUE) %>%
      visNodes(
        font = list(size = 12),
        scaling = list(min = 10, max = max_node_size),
        margin = 15
      ) %>%
      visEdges(smooth = FALSE) %>%
      visInteraction(
        zoomView = TRUE,
        dragView = TRUE,
        navigationButtons = TRUE
      ) %>%
      visEvents(selectNode = "function(params) {
    Shiny.setInputValue('selected_cluster_node', params.nodes[0]);
  }")
  })
  
  # 다운로드 버튼 UI
  output$download_cluster_network_ui <- renderUI({
    if (is.null(levene_graph_data())) {
      div(
        style = "margin-top: 20px; color: #7f8c8d; font-style: italic;",
        "No graph data available."
      )
    } else {
      downloadButton(
        "download_cluster_network",
        "Download Module Information",
        class = "btn btn-primary btn-sm",
        style = "width: 100%;"
      )
    }
  })
  
  # 다운로드 핸들러
  output$download_cluster_network <- downloadHandler(
    filename = function() {
      paste("Spatial_Analysis_Results_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      req(levene_graph_data(), fwer_results())
      
      withProgress(message = "Generating Spatial Analysis Results...", value = 0, {
        # 임시 디렉토리 생성
        temp_dir <- tempdir()
        
        # 그래프 데이터 가져오기
        graph_data <- levene_graph_data()
        g <- graph_data$cluster_graph
        cluster_sizes <- graph_data$cluster_sizes
        cluster_labels <- graph_data$clusters
        levene_genes <- graph_data$levene_genes
        
        # 데이터 유효성 검사
        if (is.null(cluster_labels) || length(cluster_labels) == 0 || any(is.na(names(cluster_labels))) || any(names(cluster_labels) == "")) {
          showNotification("Error: Invalid or empty module labels detected.", type = "error")
          return(NULL)
        }
        if (is.null(levene_genes) || length(levene_genes) == 0) {
          showNotification("Error: No valid DVGs available.", type = "error")
          return(NULL)
        }
        
        # fwer_results에서 Levene 테스트 결과 가져오기
        fwer_res <- fwer_results()
        levene_full <- fwer_res$levene_full
        levene_unique <- fwer_res$levene_unique$Gene
        
        if (is.null(levene_full) || nrow(levene_full) == 0) {
          showNotification("Error: DVG test results are empty.", type = "error")
          return(NULL)
        }
        
        # 진행바 초기화
        unique_clusters <- unique(cluster_labels)
        n_clusters <- length(unique_clusters)
        if (n_clusters == 0) {
          showNotification("Error: No valid modules found.", type = "error")
          return(NULL)
        }
        total_steps <- n_clusters + 2 # 네트워크 이미지(1) + 클러스터별 엑셀(n_clusters) + ZIP 생성(1)
        
        incProgress(1/total_steps, detail = "Generating network image...")
        
        # 최대 노드 크기와 엣지 굵기 설정
        max_node_size <- 50
        max_edge_width <- 5
        
        # 노드 크기 스케일링
        max_cluster_size <- max(cluster_sizes)
        scaled_sizes <- if (max_cluster_size > 0) {
          (cluster_sizes / max_cluster_size) * max_node_size
        } else {
          rep(max_node_size, length(cluster_sizes))
        }
        
        # 노드 데이터프레임
        nodes <- data.frame(
          id = V(g)$name,
          label = paste("Module", V(g)$name),
          size = scaled_sizes,
          color = hcl.colors(length(V(g)), palette = "Viridis"),
          title = paste("Module:", V(g)$name, "<br>Genes:", cluster_sizes)
        )
        
        # 엣지 데이터프레임
        max_edge_weight <- max(E(g)$weight, 1)
        edges <- as_data_frame(g, what = "edges") %>%
          mutate(
            width = (E(g)$weight / max_edge_weight) * max_edge_width,
            color = list(opacity = 0.4),
            title = sprintf("Edge Weight: %.3f", E(g)$weight)
          )
        
        # visNetwork 객체 생성
        network <- visNetwork(nodes, edges, height = "500px", width = "100%") %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visPhysics(enabled = TRUE) %>%
          visLayout(randomSeed = 123, improvedLayout = TRUE) %>%
          visNodes(
            font = list(size = 12),
            scaling = list(min = 10, max = max_node_size),
            margin = 15
          ) %>%
          visEdges(smooth = FALSE) %>%
          visInteraction(
            zoomView = TRUE,
            dragView = TRUE,
            navigationButtons = TRUE
          )
        
        # 네트워크 이미지 저장
        date_str <- format(Sys.Date(), "%Y%m%d")
        image_file <- paste0("cluster_network_", date_str, ".png")
        temp_html <- file.path(temp_dir, "network.html")
        temp_png <- file.path(temp_dir, image_file)
        
        visSave(network, temp_html)
        webshot(temp_html, temp_png, vwidth = 800, vheight = 500)
        
        # 클러스터별 유전자 정보 엑셀 파일 생성
        excel_files <- character(0)
        
        for (i in seq_along(unique_clusters)) {
          cluster_id <- unique_clusters[i]
          incProgress(1/total_steps, detail = sprintf("Generating Excel for Module %s...", cluster_id))
          
          # 클러스터에 속한 유전자 추출
          cluster_genes <- names(cluster_labels)[cluster_labels == cluster_id]
          
          # 유전자 정보 데이터프레임 생성
          if (length(cluster_genes) == 0) {
            showNotification(sprintf("Warning: No genes found for Module %s. Skipping Excel file creation.", cluster_id), type = "warning")
            next
          }
          
          gene_data <- data.frame(
            Gene = cluster_genes,
            P_value = NA_real_,
            Adjusted_P_value = NA_real_,
            DVG_only = FALSE
          )
          
          # Levene 테스트 결과 매핑
          for (j in seq_along(cluster_genes)) {
            gene <- cluster_genes[j]
            if (gene %in% levene_full$Gene && !is.na(gene) && gene != "") {
              gene_data$P_value[j] <- levene_full$P_value[levene_full$Gene == gene]
              gene_data$Adjusted_P_value[j] <- levene_full$`Adjusted P-value`[levene_full$Gene == gene]
              gene_data$DVG_only[j] <- gene %in% levene_unique
            }
          }
          
          # 데이터 정렬 (P_value 오름차순)
          gene_data <- gene_data %>%
            arrange(P_value)
          
          # 엑셀 파일 생성
          excel_file <- file.path(temp_dir, paste0("Module_", cluster_id, ".xlsx"))
          wb <- createWorkbook()
          addWorksheet(wb, "Genes")
          writeData(wb, "Genes", gene_data, startRow = 1, startCol = 1)
          
          # 스타일 정의
          header_style <- createStyle(
            fontColour = "#ffffff",
            fgFill = "#2980b9",
            fontSize = 11,
            fontName = "Arial",
            textDecoration = "bold",
            halign = "center",
            valign = "center"
          )
          even_row_style <- createStyle(
            fgFill = "#f8f9fa",
            halign = "center",
            valign = "center"
          )
          odd_row_style <- createStyle(
            fgFill = "#ffffff",
            halign = "center",
            valign = "center"
          )
          
          # 스타일 적용
          addStyle(wb, "Genes", style = header_style, rows = 1, cols = 1:4, gridExpand = TRUE)
          for (k in 2:(nrow(gene_data) + 1)) {
            style <- if (k %% 2 == 0) even_row_style else odd_row_style
            addStyle(wb, "Genes", style = style, rows = k, cols = 1:4, gridExpand = TRUE)
          }
          
          # 열 너비 조정
          setColWidths(wb, "Genes", cols = 1:4, widths = "auto")
          
          # 엑셀 파일 저장
          saveWorkbook(wb, excel_file, overwrite = TRUE)
          excel_files <- c(excel_files, excel_file)
        }
        
        incProgress(1/total_steps, detail = "Creating ZIP file...")
        
        # 모든 파일을 ZIP으로 압축
        files_to_zip <- c(temp_png, excel_files)
        if (length(files_to_zip) == 0) {
          showNotification("Error: No files generated for ZIP.", type = "error")
          return(NULL)
        }
        
        zip::zipr(file, files = files_to_zip, include_directories = FALSE)
        
        # 임시 파일 삭제
        unlink(c(temp_html, temp_png, excel_files))
      })
    }
  )
  
  output$cluster_genes_cards_ui <- renderUI({
    req(data_matched(), levene_graph_data(), selected_cluster_node(), fwer_results())
    graph_data <- levene_graph_data()
    cluster_id <- selected_cluster_node()
    
    if (is.null(cluster_id) || is.null(graph_data)) {
      return(div(
        style = "color: #7f8c8d; font-style: italic; text-align: center; padding: 20px;",
        "Select a module to view genes."
      ))
    }
    
    cluster_genes <- names(graph_data$clusters)[graph_data$clusters == as.numeric(cluster_id)]
    levene_fwer <- fwer_results()$levene_fwer
    levene_unique <- fwer_results()$levene_unique
    
    df <- data.frame(Gene = cluster_genes) %>%
      left_join(select(levene_fwer, Gene, Levene_Rank = `P-value rank`), by = "Gene") %>%
      left_join(select(levene_unique, Gene, Levene_Unique_Rank = `P-value rank`), by = "Gene") %>%
      mutate(
        Levene_Rank = ifelse(is.na(Levene_Rank), "Not Significant", as.character(Levene_Rank)),
        Levene_Unique_Rank = ifelse(is.na(Levene_Unique_Rank), "X", as.character(Levene_Unique_Rank))
      )
    
    search_input <- tags$input(
      type = "text",
      class = "search-bar",
      placeholder = "Search Genes...",
      oninput = "filterCards(this.value, 'cluster-genes')"
    )
    
    cards <- lapply(1:nrow(df), function(i) {
      gene <- df$Gene[i]
      levene_rank <- df$Levene_Rank[i]
      unique_rank <- df$Levene_Unique_Rank[i]
      div(
        class = "gene-card",
        `data-gene` = tolower(gene),
        onclick = sprintf("Shiny.setInputValue('selected_levene_gene_from_table', '%s')", gene),
        div(class = "gene-card-header", gene),
        div(
          class = "gene-card-content",
          tags$span(sprintf("DVG Rank: %s", levene_rank)),
          tags$span(sprintf("DVG-only Rank: %s", unique_rank))
        )
      )
    })
    
    tagList(
      search_input,
      div(class = "cards-container", id = "cluster-genes", cards)
    )
  })
  
  
  # 클러스터 노드 선택
  
  observeEvent(input$selected_cluster_node, {
    selected_cluster_node(input$selected_cluster_node)
  })
  
  # Select by id 드롭다운에서 선택했을 때도 반응하도록 추가
  observeEvent(input$cluster_network_plot_selected, {
    if (!is.null(input$cluster_network_plot_selected) && input$cluster_network_plot_selected != "") {
      selected_cluster_node(input$cluster_network_plot_selected)
    }
  })
  
  # 유전자 선택 및 분석
  observeEvent(list(input$selected_levene_gene_from_table), {
    gene <- input$selected_levene_gene_from_table
    if (is.null(gene)) return()
    
    selected_levene_gene(gene)
    
    withProgress(message = "Analyzing Selected Gene...", value = 0, {
      req(expr_data(), spots_data(), gene)
      expr <- expr_data()
      spots <- spots_data()
      
      clusters <- clusters_full()
      selected_spots <- clusters != 0
      expr_selected <- expr[, selected_spots, drop = FALSE]
      clusters_selected <- factor(clusters[selected_spots])
      
      if (!gene %in% rownames(expr_selected)) {
        showNotification(paste("Gene", gene, "not found in expression data."), type = "warning")
        return()
      }
      
      incProgress(0.5, detail = "Calculating statistics...")
      gene_data <- expr_selected[gene, ]
      mean_vals <- tapply(gene_data, clusters_selected, mean, na.rm = TRUE)
      var_vals <- tapply(gene_data, clusters_selected, var, na.rm = TRUE)
      df <- data.frame(
        Cluster = names(mean_vals),
        Mean = as.numeric(mean_vals),
        Variance = as.numeric(var_vals),
        Gene = gene
      )
      
      graph_data <- levene_graph_data()
      if (!is.null(graph_data) && gene %in% names(graph_data$clusters)) {
        df$Leiden_Cluster <- graph_data$clusters[gene]
      } else {
        df$Leiden_Cluster <- NA
      }
      
      selected_levene_gene_stats(df)
      incProgress(1.0, detail = "Complete")
    })
  })
  
  # 유전자 정보 출력
  output$levene_gene_analysis <- renderUI({
    req(data_matched())
    gene <- selected_levene_gene()
    if (is.null(gene)) {
      return(div(
        class = "no-data",
        style = "text-align: center; color: #7f8c8d; padding: 20px;",
        "No gene selected. Click a gene in the table to view details."
      ))
    }
    
    gene_stats <- selected_levene_gene_stats()
    if (is.null(gene_stats)) {
      return(div(
        class = "no-data",
        style = "text-align: center; color: #e74c3c; padding: 20px;",
        paste("No data available for gene:", gene)
      ))
    }
    
    fwer_res <- fwer_results()
    p_value <- NA
    adj_p_value <- NA
    if (!is.null(fwer_res$levene_full)) {
      levene_full <- fwer_res$levene_full
      if (gene %in% levene_full$Gene) {
        p_value <- levene_full$P_value[levene_full$Gene == gene]
        adj_p_value <- levene_full$`Adjusted P-value`[levene_full$Gene == gene]
      }
    }
    
    p_value_text <- format_scientific_notation(p_value)
    adj_p_value_text <- format_scientific_notation(adj_p_value)
    
    tagList(
      div(
        class = "gene-analysis-container",
        h3(
          tags$span("Gene for Analysis: "),
          tags$a(
            href = sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene),
            target = "_blank",
            class = "clickable-gene",
            gene
          ),
          style = "font-size: 18px; color: #2c3e50; margin-bottom: 5px;"
        ),
        p(paste("Brown-Forsythe p-value:", p_value_text), style = "font-size: 14px; color: #2c3e50; margin-bottom: 2px;"),
        p(paste("Adjusted p-value:", adj_p_value_text), style = "font-size: 14px; color: #2c3e50; margin-bottom: 15px;"),
        fluidRow(
          column(
            6, 
            plotOutput("levene_spatial_cluster_plot", height = "250px"),
            style = "padding: 10px; margin-bottom: 0;"
          ),
          column(
            6, 
            plotOutput("levene_expr_plot", height = "250px"),
            style = "padding: 10px; margin-bottom: 0;"
          )
        ),
        fluidRow(
          column(
            6, 
            uiOutput("levene_cluster_legend"),
            style = "text-align: left; padding: 10px;"
          ),
          column(
            6, 
            uiOutput("levene_expr_legend"),
            style = "text-align: right; padding: 10px;"
          )
        ),
        div(
          style = "margin: 0; padding: 0; max-height: 300px; overflow-y: auto; margin-top: 20px;",
          DTOutput("levene_stats_table", width = "100%")
        )
      )
    )
  })
  
  # 클러스터 플롯
  output$levene_spatial_cluster_plot <- renderPlot({
    req(expr_data(), coords_data(), selected_levene_gene(), selected_levene_gene_stats())
    clusters <- clusters_full()
    create_cluster_plot(coords_data(), clusters, input$spatial_point_size, input$spatial_point_alpha)
  }, height = 250)
  
  # 표현량 플롯
  output$levene_expr_plot <- renderPlot({
    req(expr_data(), coords_data(), selected_levene_gene(), selected_levene_gene_stats())
    gene <- selected_levene_gene()
    expr <- expr_data()
    clusters <- clusters_full()
    create_expression_plot(coords_data(), clusters, as.numeric(expr[gene, ]), gene,
                           input$spatial_point_size, input$spatial_point_alpha)
  }, height = 250)
  
  # 클러스터 범례
  output$levene_cluster_legend <- renderUI({
    req(expr_data(), coords_data(), selected_levene_gene(), selected_levene_gene_stats())
    create_cluster_legend(clusters_full())
  })
  
  # 표현량 범례
  output$levene_expr_legend <- renderUI({
    req(expr_data(), coords_data(), selected_levene_gene(), selected_levene_gene_stats())
    gene <- selected_levene_gene()
    expr <- expr_data()
    clusters <- clusters_full()
    df <- data.frame(Expression = as.numeric(expr[gene, ]), Cluster = clusters)
    expr_range <- range(df$Expression[clusters != "0"], na.rm = TRUE)
    create_expression_legend(expr_range)
  })
  
  output$levene_stats_table <- renderDT({
    req(selected_levene_gene_stats())
    df <- selected_levene_gene_stats() %>% 
      select(Cluster, Mean, Variance) %>%
      mutate(
        Mean = format(round(Mean, 5), nsmall = 5),
        Variance = format(round(Variance, 5), nsmall = 5)
      ) %>%
      rename(Domain = Cluster)
    datatable(
      df,
      options = list(
        dom = 't',
        paging = FALSE,
        scrollY = "200px",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames = FALSE,
      class = "compact"
    )
  })
  
  output$run_spatial_analysis_ui <- renderUI({
    is_data_matched <- data_matched()
    files_present <- !is.null(expr_data()) && !is.null(coords_data())
    
    btn <- actionButton(
      "run_spatial_analysis",
      tags$span(icon("play"), " Run Clustering"),
      class = "btn-primary btn-sm",
      width = "100%"
    )
    
    if (!files_present || !is_data_matched) {
      btn <- tagAppendAttributes(btn, disabled = "disabled")
    }
    
    btn
  })
  
  output$cluster_table_title <- renderUI({
    cluster_id <- selected_cluster_node()
    title_text <- if (is.null(cluster_id)) {
      "Genes in Selected Module"
    } else {
      paste("Genes in Module", cluster_id)
    }
    tags$h4(title_text, style = "color: #2980b9; font-size: 16px; margin-top: 20px; margin-bottom: 10px;")
  })
  
  # 유전자 검색 UI (자동완성 드롭다운)
  output$gene_search_ui <- renderUI({
    graph_data <- levene_graph_data()
    
    if (is.null(graph_data)) {
      return(NULL)
    }
    
    gene_names <- names(graph_data$clusters)
    
    tagList(
      tags$h4("Find Gene's Module", style = "color: #2980b9; font-size: 16px; margin-top: 10px; margin-bottom: 10px;"),
      selectizeInput(
        "gene_search_input",
        label = NULL,
        choices = NULL,
        selected = NULL,
        multiple = FALSE,
        width = "100%",
        options = list(
          placeholder = "Type gene name...",
          maxOptions = 10,
          create = FALSE,
          onInitialize = I('function() { this.clear(); }')
        )
      ),
      uiOutput("gene_search_result")
    )
  })
  
  observe({
    graph_data <- levene_graph_data()
    
    if (!is.null(graph_data)) {
      gene_names <- names(graph_data$clusters)
      
      updateSelectizeInput(
        session,
        "gene_search_input",
        choices = gene_names,
        selected = "",
        server = TRUE,
        options = list(
          placeholder = "Type gene name...",
          maxOptions = 10,
          create = FALSE
        )
      )
    }
  })
  
  # 유전자 검색 결과
  output$gene_search_result <- renderUI({
    graph_data <- levene_graph_data()
    search_term <- input$gene_search_input
    
    if (is.null(graph_data) || is.null(search_term) || search_term == "") {
      return(NULL)
    }
    
    cluster_labels <- graph_data$clusters
    gene_names <- names(cluster_labels)
    
    # 선택된 유전자가 클러스터에 있는지 확인
    if (search_term %in% gene_names) {
      gene <- search_term
      cluster_id <- cluster_labels[gene]
      
      # 해당 클러스터의 색상 가져오기
      unique_clusters <- unique(cluster_labels)
      n_clusters <- length(unique_clusters)
      cluster_colors <- hcl.colors(n_clusters, palette = "Viridis")
      names(cluster_colors) <- as.character(unique_clusters)
      cluster_color <- cluster_colors[as.character(cluster_id)]
      
      return(div(
        style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; border-left: 4px solid #28a745;",
        tags$span(
          style = "font-weight: bold; color: #155724;",
          gene
        ),
        tags$span(
          style = "color: #155724;",
          " is in "
        ),
        tags$span(
          style = paste0("font-weight: bold; color: ", cluster_color, "; cursor: pointer;"),
          onclick = sprintf("Shiny.setInputValue('cluster_network_plot_selected', '%s', {priority: 'event'});", cluster_id),
          paste("Module", cluster_id)
        ),
        tags$span(
          style = "color: #6c757d; font-size: 11px; margin-left: 5px;",
          "(click to view genes)"
        )
      ))
    } else {
      return(div(
        style = "background-color: #f8d7da; padding: 10px; border-radius: 5px; margin-bottom: 10px; border-left: 4px solid #dc3545;",
        tags$span(
          style = "color: #721c24;",
          paste0("Gene '", search_term, "' not found in any module.")
        )
      ))
    }
  })
  
  # Spatial Analysis 다운로드 버튼 UI
  output$download_spatial_analysis_ui <- renderUI({
    if (is.null(selected_levene_gene()) || is.null(selected_levene_gene_stats()) || 
        is.null(expr_data()) || is.null(coords_data())) {
      div(
        style = "margin-top: 20px; color: #7f8c8d; font-style: italic;",
        "No gene analysis data available."
      )
    } else {
      downloadButton(
        "download_spatial_analysis",
        "Download Gene Analysis (Image + Stats)",
        class = "btn btn-primary btn-sm",
        style = "margin-top: 20px; width: 100%;"
      )
    }
  })
  
  # Spatial Analysis 다운로드 핸들러
  output$download_spatial_analysis <- downloadHandler(
    filename = function() {
      paste("spatial_analysis_", selected_levene_gene(), "_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      if (is.null(selected_levene_gene()) || is.null(selected_levene_gene_stats()) || 
          is.null(expr_data()) || is.null(coords_data())) {
        showNotification("No gene analysis data available to download.", type = "error")
        return(NULL)
      }
      
      gene <- selected_levene_gene()
      expr <- expr_data()
      coords <- coords_data()
      stats <- selected_levene_gene_stats()
      
      fwer_res <- fwer_results()
      p_value <- NA
      adj_p_value <- NA
      
      if (!is.null(fwer_res) && !is.null(fwer_res$levene_full)) {
        levene_full <- fwer_res$levene_full
        if (gene %in% levene_full$Gene) {
          p_value <- levene_full$P_value[levene_full$Gene == gene]
          adj_p_value <- levene_full$`Adjusted P-value`[levene_full$Gene == gene]
        }
      }
      
      # 과학적 표기법 형식화 함수 (상첨자 표현식용 리스트 반환)
      format_scientific <- function(value) {
        if (length(value) == 0 || is.na(value)) return(NULL)
        sci <- format(value, digits = 5, scientific = TRUE)
        parts <- unlist(strsplit(sci, "e"))
        mantissa <- parts[1]
        exponent <- as.numeric(parts[2])
        list(mantissa = mantissa, exponent = exponent)
      }
      
      p_value_formatted <- format_scientific(p_value)
      adj_p_value_formatted <- format_scientific(adj_p_value)
      
      # 클러스터 데이터
      clusters <- clusters_full()
      cluster_levels <- levels(clusters)
      non_zero_levels <- cluster_levels[cluster_levels != "0"]
      n_colors <- length(non_zero_levels)
      colors <- hcl.colors(n_colors, palette = "Set2")
      color_mapping <- setNames(c("#bdc3c7", colors), c("0", non_zero_levels))
      
      # 플롯 1: 클러스터 분포 (Spatial Analysis 탭의 플롯 사용)
      p1 <- ggplot(coords, aes(x = X, y = Y)) +
        geom_point(data = coords[clusters == "0", ], 
                   color = "#bdc3c7", shape = 1, size = 5 * input$spatial_point_size, alpha = input$spatial_point_alpha) +
        geom_point(data = coords[clusters != "0", ], 
                   aes(color = factor(clusters[clusters != "0"])), 
                   shape = 16, size = 5 * input$spatial_point_size, alpha = input$spatial_point_alpha) +
        scale_color_manual(values = color_mapping) +
        labs(title = "Tissue Domains", x = "X Coordinate", y = "Y Coordinate") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none",
          panel.grid = element_blank()
        ) +
        coord_fixed()
      
      # 플롯 2: 유전자 발현 (Spatial Analysis 탭의 플롯 사용)
      p2 <- ggplot(coords, aes(x = X, y = Y)) +
        geom_point(data = coords[clusters == "0", ], 
                   color = "#bdc3c7", shape = 1, size = 5 * input$spatial_point_size, alpha = input$spatial_point_alpha) +
        geom_point(data = coords[clusters != "0", ], 
                   aes(color = expr[gene, clusters != "0"]), 
                   shape = 16, size = 5 * input$spatial_point_size, alpha = input$spatial_point_alpha) +
        scale_color_viridis_c(option = "viridis", name = "Expression") +
        labs(title = paste("Expression of", gene), x = "X Coordinate", y = "Y Coordinate") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none",
          panel.grid = element_blank()
        ) +
        coord_fixed()
      
      # 플롯 결합
      plot_row <- plot_grid(p1, p2, ncol = 2, align = "hv", axis = "tblr")
      plot_height <- 300  # 플롯 높이 (포인트 단위)
      plot_width <- 600   # 플롯 너비 (포인트 단위)
      
      # 클러스터 범례
      items <- lapply(1:n_colors, function(i) {
        list(color = colors[i], label = non_zero_levels[i])
      })
      rows <- split(items, ceiling(seq_along(items) / 5))
      legend_height_cluster <- length(rows) * 30  # 각 행당 30pt
      legend_width_cluster <- 300  # 클러스터 범례 너비
      
      png(tempfile())
      pushViewport(viewport(layout = grid.layout(nrow = length(rows), ncol = 5, widths = unit(0.2, "npc"))))
      for (i in seq_along(items)) {
        row <- ceiling(i / 5)
        col <- (i - 1) %% 5 + 1
        pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
        grid.circle(x = 0.2, y = 0.5, r = 0.1, gp = gpar(fill = items[[i]]$color, col = NA))
        grid.text(label = items[[i]]$label, x = 0.4, y = 0.5, just = "left", gp = gpar(fontsize = 12))
        popViewport()
      }
      legend_grob_cluster <- grid.grab()
      dev.off()
      
      # 발현 범례
      expr_range <- range(expr[gene, clusters != "0"], na.rm = TRUE)
      min_val <- expr_range[1]
      max_val <- expr_range[2]
      legend_height_expr <- 50  # 발현 범례 높이
      legend_width_expr <- 150  # 발현 범례 너비
      
      png(tempfile())
      gradient_colors <- viridis::viridis(256, option = "viridis")
      gradient_matrix <- matrix(gradient_colors, nrow = 1, ncol = 256)
      pushViewport(viewport(width = unit(150, "points"), height = unit(12, "points")))
      grid.raster(gradient_matrix, width = 1, height = 1, interpolate = TRUE)
      gradient_grob <- grid.grab()
      dev.off()
      
      png(tempfile())
      pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, 
                                                 heights = unit.c(unit(20, "points"), unit(30, "points")))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.text("Expression Level", x = 0, y = 0.5, just = "left", gp = gpar(fontsize = 12, col = "#2c3e50"))
      popViewport()
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, 
                                                 heights = unit.c(unit(12, "points"), unit(18, "points")))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.draw(gradient_grob)
      popViewport()
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      grid.text(sprintf("%.2f", min_val), x = 0, y = 0.5, just = "left", gp = gpar(fontsize = 12, col = "#2c3e50"))
      grid.text(sprintf("%.2f", max_val), x = 1, y = 0.5, just = "right", gp = gpar(fontsize = 12, col = "#2c3e50"))
      popViewport()
      popViewport()
      expr_legend_grob <- grid.grab()
      dev.off()
      
      # 텍스트 영역 크기 계산
      text_height <- 100  # 텍스트 영역 높이 (제목 및 p-value 등)
      text_width <- 600   # 텍스트 영역 너비
      
      # 전체 크기 계산 (표 제외)
      total_width <- max(plot_width, text_width, legend_width_cluster + legend_width_expr + 20)
      legend_height <- max(legend_height_cluster, legend_height_expr) + 20  # 범례 높이 + 여유
      total_height <- text_height + plot_height + legend_height + 30  # 표 제외, 여유 공간 조정
      
      # 사용자 지정 파일 이름 생성
      date_str <- format(Sys.Date(), "%Y%m%d")
      image_file <- paste0("spatial_analysis_image_", gene, "_", date_str, ".png")
      excel_file <- paste0("spatial_analysis_stats_", gene, "_", date_str, ".xlsx")
      
      # 임시 파일 경로 생성
      image_path <- file.path(tempdir(), image_file)
      excel_path <- file.path(tempdir(), excel_file)
      
      # 이미지 생성
      CairoPNG(image_path, width = total_width, height = total_height, units = "pt", dpi = 300)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, 
                                                 heights = unit.c(
                                                   unit(text_height, "points"),
                                                   unit(plot_height, "points"),
                                                   unit(legend_height, "points")
                                                 ))))
      
      # 1. 텍스트 영역
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      grid.text(paste("Gene for Analysis:", gene), x = 0.1, y = 0.9, just = "left", 
                gp = gpar(fontsize = 18, fontface = "bold", col = "#2c3e50"))
      if (!is.null(p_value_formatted)) {
        grid.text(label = bquote("Brown-Forsythe p-value: " ~ .(p_value_formatted$mantissa) ~ " x " ~ 10^.(p_value_formatted$exponent)), 
                  x = 0.1, y = 0.6, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      } else {
        grid.text("Brown-Forsythe p-value: NA", x = 0.1, y = 0.6, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      }
      if (!is.null(adj_p_value_formatted)) {
        grid.text(label = bquote("Adjusted p-value: " ~ .(adj_p_value_formatted$mantissa) ~ " x " ~ 10^.(adj_p_value_formatted$exponent)), 
                  x = 0.1, y = 0.3, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      } else {
        grid.text("Adjusted p-value: NA", x = 0.1, y = 0.3, just = "left", 
                  gp = gpar(fontsize = 14, col = "#2c3e50"))
      }
      popViewport()
      
      # 2. 플롯 영역
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      grid.draw(ggplotGrob(plot_row))
      popViewport()
      
      # 3. 범례 영역 (좌우 배치)
      pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
      pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, 
                                                 widths = unit.c(unit(legend_width_cluster, "points"), unit(legend_width_expr, "points")))))
      pushViewport(viewport(layout.pos.col = 1))
      grid.draw(legend_grob_cluster)
      popViewport()
      pushViewport(viewport(layout.pos.col = 2))
      grid.draw(expr_legend_grob)
      popViewport()
      popViewport()
      
      dev.off()
      
      # 엑셀 파일 생성 및 스타일링
      stats_table <- stats %>%
        select(Cluster, Mean, Variance) %>% # Leiden_Cluster 열 제외
        mutate(
          Mean = format(round(Mean, 5), nsmall = 5),
          Variance = format(round(Variance, 5), nsmall = 5)
        ) %>%
        rename(Domain = Cluster)
      wb <- createWorkbook()
      addWorksheet(wb, "Stats")
      writeData(wb, "Stats", stats_table, startRow = 1, startCol = 1)
      
      # 스타일 정의
      header_style <- createStyle(
        fontColour = "#ffffff",
        fgFill = "#2980b9",
        fontSize = 11,
        fontName = "Arial",
        textDecoration = "bold",
        halign = "center",
        valign = "center"
      )
      even_row_style <- createStyle(
        fgFill = "#f8f9fa",
        halign = "center",
        valign = "center"
      )
      odd_row_style <- createStyle(
        fgFill = "#ffffff",
        halign = "center",
        valign = "center"
      )
      
      # 스타일 적용
      addStyle(wb, "Stats", style = header_style, rows = 1, cols = 1:3, gridExpand = TRUE)
      for (i in 2:(nrow(stats_table) + 1)) {
        style <- if (i %% 2 == 0) even_row_style else odd_row_style
        addStyle(wb, "Stats", style = style, rows = i, cols = 1:3, gridExpand = TRUE)
      }
      
      # 열 너비 조정
      setColWidths(wb, "Stats", cols = 1:3, widths = "auto")
      
      # 엑셀 파일 저장
      saveWorkbook(wb, excel_path, overwrite = TRUE)
      
      # ZIP 파일 생성
      zip::zipr(file, files = c(image_path, excel_path), include_directories = FALSE)
      
      # 임시 파일 삭제
      unlink(image_path)
      unlink(excel_path)
    }
  )
  ##############################################################################
  
  # 클러스터 다운로드 핸들러
  output$download_clusters <- downloadHandler(
    filename = function() {
      date_str <- format(Sys.Date(), "%Y%m%d")
      paste0("domain_assignments_", date_str, ".csv")
    },
    content = function(file) {
      # 필요한 데이터 가져오기
      req(expr_data(), coords_data(), clusters_full())
      
      coords <- coords_data()
      clusters <- clusters_full()
      
      # spot 이름 가져오기 (expression matrix의 column 순서와 동일)
      spot_names <- colnames(expr_data())
      
      # 클러스터 정보를 데이터프레임으로 생성
      # spot 순서는 expression matrix의 column 순서와 일치
      cluster_df <- data.frame(
        Spot = spot_names,
        Domain = as.character(clusters),
        stringsAsFactors = FALSE
      )
      
      # CSV 파일로 저장
      write.csv(cluster_df, file, row.names = FALSE)
    }
  )
  
}

# Gene 이름 클릭시 빈탭 열리는거 해결하기 위해 외부 뷰어에서 열기
options(shiny.launch.browser = TRUE)
# Shiny 앱 실행
shinyApp(ui = ui, server = server)
