using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace GammaTextScope
{
    public static class ChartPrinter
    {
        public static void PrintChart(List<KeyValuePair<string, double>> series, string chartTitle, string xAxisLabel, string yAxisLabel, string filePath)
        {
            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial", PlotType = PlotType.XY };
            model.Axes.Add(new LinearAxis { Title = xAxisLabel, Position = AxisPosition.Bottom });
            model.Axes.Add(new LinearAxis { Title = yAxisLabel, Position = AxisPosition.Left });

            var groupings = DistributionHelper.CountNumberGroupings(series.Select(x => x.Value).OrderBy(x => x).ToList(), 0.1);
            var lineSeries = new OxyPlot.Series.AreaSeries();

            foreach (var group in groupings)
                lineSeries.Points.Add(new DataPoint(group.Key, group.Value));

            model.Series.Add(lineSeries);

            var exporter = new PdfExporter { Width = 1000, Height = 500 };
            using (var stream = File.Create(@filePath))
            {
                exporter.Export(model, stream);
            }
        }

        public static void PrintChartUnsorted(List<KeyValuePair<double, double>> series, string chartTitle, string xAxisLabel, string yAxisLabel, string filePath, double abosouteMinimum, List<string> xAxixSlabels = null, double majorStep = 2)
        {
            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial", PlotType = PlotType.XY };
            model.Axes.Add(new LinearAxis { Title = yAxisLabel, Position = AxisPosition.Left, MajorStep = majorStep });

            if (xAxixSlabels != null)
            {
                var xAxis = new CategoryAxis { Title = xAxisLabel, Position = AxisPosition.Bottom, AbsoluteMinimum = abosouteMinimum };
                if (xAxixSlabels != null)
                    foreach (var label in xAxixSlabels)
                        xAxis.Labels.Add(label);

                model.Axes.Add(xAxis);
            }


            var lineSeries = new OxyPlot.Series.AreaSeries();

            foreach (var group in series)
                lineSeries.Points.Add(new DataPoint(group.Key, group.Value));

            model.Series.Add(lineSeries);

            var exporter = new PdfExporter { Width = 1000, Height = 500 };
            using (var stream = File.Create(@filePath))
            {
                exporter.Export(model, stream);
            }
        }

        public static void PrintChartBig(List<KeyValuePair<double, double>> series, List<string> xAxixSlabels, string chartTitle, string xAxisLabel, string yAxisLabel, string filePath)
        {
            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial", PlotType = PlotType.XY };
            model.Axes.Add(new LinearAxis { Title = yAxisLabel, Position = AxisPosition.Left, MajorStep = 2 });

            var xAxis = new CategoryAxis { Title = xAxisLabel, Position = AxisPosition.Bottom, AbsoluteMinimum = 1 };
            foreach (var label in xAxixSlabels)
                xAxis.Labels.Add(label);

            model.Axes.Add(xAxis);

            var lineSeries = new OxyPlot.Series.AreaSeries();

            foreach (var group in series)
                lineSeries.Points.Add(new DataPoint(group.Key, group.Value));

            model.Series.Add(lineSeries);

            var exporter = new PdfExporter { Width = 2500, Height = 500 };
            using (var stream = File.Create(@filePath))
            {
                exporter.Export(model, stream);
            }
        }

        public static void PrintHorizontalBarChart(List<double> values, List<string> labels, string chartTitle, string valueLabel, string categoryLabel, string filePath)
        {
            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial" };

            // 1. Category Axis (Y-Axis) - Chromosomes
            // We reverse the order so "Chr1" appears at the TOP, matching the table.
            var categoryAxis = new CategoryAxis
            {
                Position = AxisPosition.Left,
                Title = categoryLabel,
                StartPosition = 1,
                EndPosition = 0 // Flips axis so index 0 is at the top
            };

            // Add labels in the original order (the Axis flip handles the visual position)
            foreach (var label in labels)
                categoryAxis.Labels.Add(label);

            model.Axes.Add(categoryAxis);

            // 2. Value Axis (X-Axis) - The Mean values
            var valueAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = valueLabel,
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot
            };
            model.Axes.Add(valueAxis);

            // 3. Create Horizontal Bars
            var barSeries = new OxyPlot.Series.BarSeries
            {
                Title = chartTitle,
                StrokeColor = OxyColors.Black,
                StrokeThickness = 1,
                FillColor = OxyColors.CornflowerBlue // Nice distinct color
            };

            foreach (var v in values)
            {
                barSeries.Items.Add(new BarItem { Value = v });
            }

            model.Series.Add(barSeries);

            // 4. Export to PDF (Taller height to fit all chromosomes cleanly)
            var exporter = new PdfExporter { Width = 800, Height = 1000 };
            using (var stream = System.IO.File.Create(filePath))
            {
                exporter.Export(model, stream);
            }
        }

        /// <summary>
        /// Grouped horizontal bar chart with multiple series sharing a common category axis.
        /// Each entry in seriesData adds a named bar series; labels map to the category axis (Y).
        /// </summary>
        public static void PrintGroupedHorizontalBarChart(
            List<(string Name, List<double> Values, OxyColor Color)> seriesData,
            List<string> labels,
            string chartTitle,
            string valueLabel,
            string categoryLabel,
            string filePath)
        {
            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial" };

            var categoryAxis = new CategoryAxis
            {
                Position = AxisPosition.Left,
                Title = categoryLabel,
                StartPosition = 1,
                EndPosition = 0
            };
            foreach (var label in labels)
                categoryAxis.Labels.Add(label);
            model.Axes.Add(categoryAxis);

            var valueAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = valueLabel,
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot
            };
            model.Axes.Add(valueAxis);

            foreach (var (name, values, color) in seriesData)
            {
                var barSeries = new OxyPlot.Series.BarSeries
                {
                    Title = name,
                    FillColor = color,
                    StrokeColor = OxyColors.Black,
                    StrokeThickness = 1
                };
                foreach (var v in values)
                    barSeries.Items.Add(new BarItem { Value = v });
                model.Series.Add(barSeries);
            }

            model.IsLegendVisible = true;
            model.LegendPosition = OxyPlot.LegendPosition.BottomCenter;

            var exporter = new PdfExporter { Width = 1000, Height = 1000 };
            using (var stream = File.Create(filePath))
            {
                exporter.Export(model, stream);
            }
        }
        /// <summary>
        /// Line chart of stop codon delta per window, with colored background bands
        /// showing classification (alfa=green, zipf=blue, ikke=gray) for Normal (top band)
        /// and Tumor (bottom band). Optional chromosome boundary lines for genome-wide view.
        /// </summary>
        public static void PrintStopCodonAlphaLineChart(
            List<long> deltas,
            List<string> normalClassifications,
            List<string> tumorClassifications,
            string chartTitle,
            string xAxisLabel,
            string yAxisLabel,
            string filePath,
            List<(int Index, string Label)> chromosomeBoundaries,
            string subtitle = null,
            List<(int StartIndex, int EndIndex)> centromereRanges = null)
        {
            bool isGenomeWide = chromosomeBoundaries != null && chromosomeBoundaries.Count > 0;
            int width = isGenomeWide ? 3000 : 1200;
            int height = 600;

            var model = new PlotModel { Title = chartTitle, DefaultFont = "Arial" };
            if (!string.IsNullOrEmpty(subtitle))
            {
                model.Subtitle = subtitle;
                model.SubtitleFontSize = 10;
                model.SubtitleColor = OxyColors.Gray;
            }

            var xAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = xAxisLabel,
                Minimum = -0.5,
                Maximum = deltas.Count - 0.5
            };
            model.Axes.Add(xAxis);

            // Find Y range for the delta line
            double yMin = deltas.Count > 0 ? deltas.Min() : -1;
            double yMax = deltas.Count > 0 ? deltas.Max() : 1;
            double yPad = Math.Max(1, (yMax - yMin) * 0.15);

            // Reserve space: both bands + decay marker rows go BELOW the delta line
            double bandHeight = (yMax - yMin) * 0.14;
            if (bandHeight < 2) bandHeight = Math.Max(2, yPad * 0.5);
            double markerRowHeight = bandHeight * 0.6;
            double gap = yPad * 0.15; // small gap between elements

            // Layout below the data (top-to-bottom):
            //   yMin - yPad  →  Tumor band top
            //   Tumor band (height = bandHeight)
            //   Tumor decay dots row
            //   gap
            //   Normal band (height = bandHeight)
            //   Normal decay dots row
            double tumorBandTop = yMin - yPad;
            double tumorBandBottom = tumorBandTop - bandHeight;
            double tumorMarkerY = tumorBandBottom - markerRowHeight * 0.5;

            double normalBandTop = tumorBandBottom - markerRowHeight - gap;
            double normalBandBottom = normalBandTop - bandHeight;
            double normalMarkerY = normalBandBottom - markerRowHeight * 0.5;

            double totalTop = yMax + yPad;
            double totalBot = normalBandBottom - markerRowHeight - yPad * 0.3;

            var yAxis = new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = yAxisLabel,
                Minimum = totalBot,
                Maximum = totalTop,
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot
            };
            model.Axes.Add(yAxis);

            // Draw classification rectangles
            for (int i = 0; i < deltas.Count; i++)
            {
                // Normal classification band (top)
                model.Annotations.Add(new RectangleAnnotation
                {
                    MinimumX = i - 0.5,
                    MaximumX = i + 0.5,
                    MinimumY = normalBandBottom,
                    MaximumY = normalBandTop,
                    Fill = ClassColor(normalClassifications[i]),
                    StrokeThickness = 0
                });

                // Tumor classification band (bottom)
                model.Annotations.Add(new RectangleAnnotation
                {
                    MinimumX = i - 0.5,
                    MaximumX = i + 0.5,
                    MinimumY = tumorBandBottom,
                    MaximumY = tumorBandTop,
                    Fill = ClassColor(tumorClassifications[i]),
                    StrokeThickness = 0
                });
            }

            // Band labels
            model.Annotations.Add(new TextAnnotation
            {
                TextPosition = new DataPoint(-1.5, (normalBandBottom + normalBandTop) / 2),
                Text = "Normal",
                StrokeThickness = 0,
                FontSize = 9
            });
            model.Annotations.Add(new TextAnnotation
            {
                TextPosition = new DataPoint(-1.5, (tumorBandBottom + tumorBandTop) / 2),
                Text = "Tumour",
                StrokeThickness = 0,
                FontSize = 9
            });

            // Qualification decay delta markers: sliding window of 10, mark where
            // Normal and Tumor differ by 3+ qualified windows
            const int slidingSize = 10;
            const int minDelta = 3;
            double dotSize = 5; // fixed size, independent of data range
            AddDecayDeltaMarkers(model, normalClassifications, tumorClassifications,
                normalMarkerY, tumorMarkerY, dotSize, slidingSize, minDelta);

            // Delta line
            var lineSeries = new OxyPlot.Series.LineSeries
            {
                Title = "Delta (Tumour − Normal)",
                Color = OxyColors.Black,
                StrokeThickness = 1.5
            };
            for (int i = 0; i < deltas.Count; i++)
                lineSeries.Points.Add(new DataPoint(i, deltas[i]));
            model.Series.Add(lineSeries);

            // Zero line
            model.Annotations.Add(new LineAnnotation
            {
                Type = LineAnnotationType.Horizontal,
                Y = 0,
                Color = OxyColors.Red,
                StrokeThickness = 0.8,
                LineStyle = LineStyle.Dash
            });

            // Centromere bands — subtle background rectangles spanning the full Y range
            // Minimum visual width ensures even single-window centromeres are visible
            if (centromereRanges != null && centromereRanges.Count > 0)
            {
                var centromereColor = OxyColor.FromArgb(40, 255, 180, 0); // very subtle amber/yellow
                double minVisualWidth = Math.Max(3, deltas.Count * 0.004); // at least 0.4% of chart width
                foreach (var (cStart, cEnd) in centromereRanges)
                {
                    double actualWidth = (cEnd + 0.5) - (cStart - 0.5);
                    double drawWidth = Math.Max(actualWidth, minVisualWidth);
                    double center = (cStart + cEnd) / 2.0;

                    model.Annotations.Add(new RectangleAnnotation
                    {
                        MinimumX = center - drawWidth / 2,
                        MaximumX = center + drawWidth / 2,
                        MinimumY = totalBot,
                        MaximumY = totalTop,
                        Fill = centromereColor,
                        StrokeThickness = 0,
                        Layer = AnnotationLayer.BelowSeries
                    });
                }
            }

            // Chromosome boundaries (genome-wide chart)
            if (isGenomeWide)
            {
                foreach (var (idx, label) in chromosomeBoundaries)
                {
                    model.Annotations.Add(new LineAnnotation
                    {
                        Type = LineAnnotationType.Vertical,
                        X = idx - 0.5,
                        Color = OxyColors.LightGray,
                        StrokeThickness = 0.5,
                        LineStyle = LineStyle.Dash
                    });
                    model.Annotations.Add(new TextAnnotation
                    {
                        TextPosition = new DataPoint(idx + 1, totalTop - bandHeight * 0.3),
                        Text = label.Replace("chr", ""),
                        StrokeThickness = 0,
                        FontSize = 7,
                        TextColor = OxyColors.Gray
                    });
                }
            }

            // Legend rectangles for classification colors
            // Add invisible scatter series to get legend entries
            model.Series.Add(new OxyPlot.Series.ScatterSeries { Title = "Alfa", MarkerFill = OxyColor.FromArgb(220, 220, 40, 40), MarkerType = MarkerType.Square, MarkerSize = 0 });
            model.Series.Add(new OxyPlot.Series.ScatterSeries { Title = "Zipf", MarkerFill = OxyColor.FromArgb(220, 60, 120, 230), MarkerType = MarkerType.Square, MarkerSize = 0 });
            model.Series.Add(new OxyPlot.Series.ScatterSeries { Title = "Unqualified", MarkerFill = OxyColor.FromArgb(150, 190, 190, 190), MarkerType = MarkerType.Square, MarkerSize = 0 });
            model.IsLegendVisible = true;
            model.LegendPosition = LegendPosition.BottomCenter;
            model.LegendOrientation = LegendOrientation.Horizontal;

            var exporter = new PdfExporter { Width = width, Height = height };
            using (var stream = File.Create(filePath))
            {
                exporter.Export(model, stream);
            }
        }

        private static OxyColor ClassColor(string classification)
        {
            switch (classification)
            {
                case "alfa": return OxyColor.FromArgb(220, 220, 40, 40);   // red
                case "zipf": return OxyColor.FromArgb(220, 60, 120, 230);  // vivid blue
                default:     return OxyColor.FromArgb(150, 190, 190, 190); // light gray
            }
        }

        /// <summary>
        /// Scans Normal vs Tumor classifications with a sliding window.
        /// Where one side has 3+ more qualified windows (alfa+zipf) than the other,
        /// places a decay marker:
        ///   - Dot above Normal band = Normal more qualified (tumor decayed)
        ///   - Dot below Tumor band  = Tumor more qualified (normal decayed)
        /// Green dot if the advantage is alfa-dominated, blue if zipf-dominated.
        /// </summary>
        private static void AddDecayDeltaMarkers(
            PlotModel model,
            List<string> normalClassifications,
            List<string> tumorClassifications,
            double normalMarkerY,
            double tumorMarkerY,
            double dotSize,
            int slidingSize,
            int minDelta)
        {
            int n = Math.Min(normalClassifications.Count, tumorClassifications.Count);
            if (n < slidingSize) return;

            for (int center = 0; center < n; center++)
            {
                int start = Math.Max(0, center - slidingSize / 2);
                int end = Math.Min(n, start + slidingSize);
                if (end - start < slidingSize && end == n) start = Math.Max(0, end - slidingSize);

                int nAlfa = 0, nZipf = 0, tAlfa = 0, tZipf = 0;
                for (int j = start; j < end; j++)
                {
                    if (normalClassifications[j] == "alfa") nAlfa++;
                    else if (normalClassifications[j] == "zipf") nZipf++;
                    if (tumorClassifications[j] == "alfa") tAlfa++;
                    else if (tumorClassifications[j] == "zipf") tZipf++;
                }

                int normalQualified = nAlfa + nZipf;
                int tumorQualified = tAlfa + tZipf;
                int delta = normalQualified - tumorQualified;

                if (Math.Abs(delta) >= minDelta)
                {
                    // Determine dominant class on the side that has more
                    bool alfaDominant;
                    double markerY;

                    if (delta > 0)
                    {
                        // Normal has more qualified → tumor decayed → marker above Normal band
                        alfaDominant = nAlfa > nZipf;
                        markerY = normalMarkerY;
                    }
                    else
                    {
                        // Tumor has more qualified → normal decayed → marker below Tumor band
                        alfaDominant = tAlfa > tZipf;
                        markerY = tumorMarkerY;
                    }

                    OxyColor dotColor = alfaDominant
                        ? OxyColor.FromArgb(230, 220, 40, 40)   // red (alfa decay)
                        : OxyColor.FromArgb(230, 60, 120, 230);  // vivid blue (zipf decay)

                    model.Annotations.Add(new PointAnnotation
                    {
                        X = center,
                        Y = markerY,
                        Fill = dotColor,
                        Size = dotSize,
                        StrokeThickness = 0,
                        Shape = MarkerType.Circle
                    });
                }
            }
        }
    }
}
