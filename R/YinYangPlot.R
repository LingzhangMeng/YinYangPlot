#' Create a Volcano Plot for Differential Expression Analysis
#'
#' Enhanced Volcano Plot for RNA-seq DE Results
#'
#' @param df Data frame with DE results
#' @param log2fc_col Column name for log2 fold change (default: "log2FoldChange")
#' @param padj_col Column name for adjusted p-value (default: "padj")
#' @param gene_col Column name for gene symbols (default: "gene_symbol")
#' @param fc_cutoff Threshold for fold change significance (default: 0.5)
#' @param p_cutoff Threshold for adjusted p-value significance (default: 0.05)
#' @param top_n_labels Number of top genes to label (default: 20)
#' @param boxed_labels Logical, whether labels have boxes (default: TRUE)
#' @param col_up Color for upregulated genes (default: "red")
#' @param col_down Color for downregulated genes (default: "blue")
#' @param label_size Font size for gene labels (default: 3)
#' @return ggplot2 object
#' @importFrom utils head
#' @export
yy_volcano <- function(
    df,
    log2fc_col = "log2FoldChange",
    padj_col   = "padj",
    gene_col   = "gene_symbol",
    fc_cutoff  = 0.5,
    p_cutoff   = 0.05,
    top_n_labels = 20,
    boxed_labels = TRUE,
    col_up     = "red",       # NEW: color for upregulated
    col_down   = "blue",      # NEW: color for downregulated
    label_size = 3             # NEW: font size for gene labels
) {
  df <- df %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        (!!as.name(padj_col) < p_cutoff & !!as.name(log2fc_col) > fc_cutoff) ~ "Up",
        (!!as.name(padj_col) < p_cutoff & !!as.name(log2fc_col) < -fc_cutoff) ~ "Down",
        TRUE ~ "NotSignificant"
      )
    )

  # Count DEGs
  n_up   <- sum(df$sig == "Up")
  n_down <- sum(df$sig == "Down")

  # Volcano plot
  p <- ggplot(df, aes(x = !!as.name(log2fc_col), y = -log10(!!as.name(padj_col)))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = col_up, "Down" = col_down, "NotSignificant" = "grey70")) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    theme_bw() +
    theme(
      panel.background = element_blank(),    # no background
      panel.grid = element_blank(),          # no grid
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "transparent"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Volcano Plot",
      subtitle = paste0(
        "DEGs: Up=", n_up, ", Down=", n_down
      ),
      color = "Significance",
      x = "log2 Fold Change",
      y = "-log10(adj.p-value)"
    )

  # Add labels
  if (top_n_labels > 0) {
    top_genes <- df %>%
      dplyr::filter(sig != "NotSignificant") %>%
      dplyr::arrange(!!as.name(padj_col)) %>%
      head(top_n_labels)

    p <- p +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = !!as.name(gene_col)),
        size = label_size,
        box.padding = 0.5,
        segment.color = "grey50",
        max.overlaps = Inf
      )
  }

  if (boxed_labels) {
    p <- p + ggrepel::geom_label_repel(
      data = top_genes,
      aes(label = !!as.name(gene_col)),
      size = label_size,
      box.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = Inf
    )
  } else {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = !!as.name(gene_col)),
      size = label_size,
      box.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = Inf
    )
  }

  return(p)
}






#' Create Yin-Yang Plot
#'
#' @description
#' Generates a circular yin-yang plot for visualizing differential gene expression.
#' In this plot, upregulated genes (Yang) appear on the right side, downregulated
#' genes (Yin) appear on the left side. Distance from the center represents
#' absolute log2 fold change, and point size represents statistical significance.
#'
#' @param data A data frame containing differential expression results.
#'   Must contain columns specified in `gene_id_col`, `gene_symbol_col`,
#'   `log2fc_col`, `padj_col`, and `significance_col`.
#' @param gene_id_col Name of column containing unique gene identifiers.
#'   Default: "gene_id"
#' @param gene_symbol_col Name of column containing gene symbols for labeling.
#'   Default: "gene_symbol"
#' @param log2fc_col Name of column containing log2 fold change values.
#'   Default: "log2FoldChange"
#' @param padj_col Name of column containing adjusted p-values.
#'   Default: "padj"
#' @param significance_col Name of column indicating significance direction.
#'   Should contain "Upregulated" or "Downregulated" for significant genes.
#'   Default: "significance"
#' @param show_all Logical. If `TRUE`, show all significant genes.
#'   If `FALSE`, show top `top_n_up` upregulated and `top_n_down` downregulated genes.
#'   Default: `TRUE`
#' @param top_n_up Number of top upregulated genes to show when `show_all = FALSE`.
#'   Default: 5
#' @param top_n_down Number of top downregulated genes to show when `show_all = FALSE`.
#'   Default: 5
#' @param gene_label_size Font size for gene labels. If `NULL`, calculated automatically
#'   based on available space and number of labels. Default: `NULL`
#' @param point_size_range Range for point sizes (minimum and maximum).
#'   Default: `c(1, 7)`
#' @param background_alpha Alpha transparency for background semicircles (0-1).
#'   Default: 0.2
#' @param show_legend Logical. Show legend on right side? Default: `TRUE`
#' @param label_top_n Number of top genes to label by significance.
#'   Default: 20
#' @param title Plot title. Default: "Yin-Yang Gene Regulation Map"
#' @param subtitle Plot subtitle. If `NULL`, generated automatically.
#'   Default: `NULL`
#' @param yin_color Background color for yin (downregulated) semicircle.
#'   Default: "#90EE90" (light green)
#' @param yang_color Background color for yang (upregulated) semicircle.
#'   Default: "#FFB6C1" (light pink)
#' @param point_yin_color Point color for downregulated genes.
#'   Default: "#006400" (dark green)
#' @param point_yang_color Point color for upregulated genes.
#'   Default: "#8B0000" (dark red)
#' @param border_color Color for plot border. If `NULL`, no border.
#'   Default: `NULL`
#' @param border_size Size of plot border (if `border_color` not `NULL`).
#'   Default: 1
#' @param min_radius Minimum radius for points (0-1). Default: 0.3
#' @param max_radius Maximum radius for points (0-1). Default: 0.85
#' @param radius_power Power transformation for radius scaling.
#'   1 = linear, 2 = quadratic, 0.5 = square root. Default: 1
#' @param plot_device_width Width of plotting device in inches for auto-sizing.
#'   If `NULL`, uses default. Default: `NULL`
#' @param plot_device_height Height of plotting device in inches for auto-sizing.
#'   If `NULL`, uses default. Default: `NULL`
#' @param seed Random seed for reproducible gene placement. Default: 42
#'
#' @return A ggplot2 object representing the yin-yang plot.
#'
#' @details
#' The yin-yang plot provides an intuitive visualization of differential gene expression:
#' \itemize{
#'   \item \strong{Yin (left side)}: Downregulated genes
#'   \item \strong{Yang (right side)}: Upregulated genes
#'   \item \strong{Distance from center}: Absolute log2 fold change (|log2FC|)
#'   \item \strong{Point size}: Statistical significance (-log10(adjusted p-value))
#'   \item \strong{Point color}: Regulation direction (green = down, red = up)
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(de_results)
#'
#' # Basic yin-yang plot
#' p <- yinyang(data = de_results)
#' print(p)
#'
#' # Customized plot
#' p_custom <- yinyang(
#'   data = de_results,
#'   show_all = FALSE,
#'   top_n_up = 10,
#'   top_n_down = 10,
#'   yin_color = "#E8F4F8",
#'   yang_color = "#FDEDEC",
#'   point_yin_color = "#2980B9",
#'   point_yang_color = "#E74C3C",
#'   title = "Differential Gene Expression",
#'   label_top_n = 15
#' )
#' print(p_custom)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @import scales
#' @importFrom stringr str_sub
#' @importFrom grDevices dev.list dev.size
#' @importFrom stats quantile
#' @importFrom stats runif quantile
#' @importFrom utils head
#' @importFrom grDevices dev.list dev.size
#' @export
yinyang <- function(
    data,
    gene_id_col = "gene_id",
    gene_symbol_col = "gene_symbol",
    log2fc_col = "log2FoldChange",
    padj_col = "padj",
    significance_col = "significance",
    show_all = TRUE,
    top_n_up = 5,
    top_n_down = 5,
    gene_label_size = NULL,
    point_size_range = c(1, 7),
    background_alpha = 0.2,
    show_legend = TRUE,
    label_top_n = 20,
    title = "Yin-Yang Gene Regulation Map",
    subtitle = NULL,
    yin_color = "#90EE90",
    yang_color = "#FFB6C1",
    point_yin_color = "#006400",
    point_yang_color = "#8B0000",
    border_color = NULL,
    border_size = 1,
    min_radius = 0.3,
    max_radius = 0.85,
    radius_power = 1,
    plot_device_width = NULL,
    plot_device_height = NULL,
    seed = 42
) {

  # Validate required packages
  required_packages <- c("ggplot2", "dplyr", "ggrepel", "scales")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Required package(s) not installed: ", paste(missing_packages, collapse = ", "))
  }

  # Validate input data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Check required columns
  required_cols <- c(gene_id_col, gene_symbol_col, log2fc_col, padj_col, significance_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Filter for significant genes only
  sig_genes <- data %>%
    dplyr::filter(!!as.name(significance_col) %in% c("Upregulated", "Downregulated"))

  if (nrow(sig_genes) == 0) {
    warning("No significant genes found. Please check your significance column.")
    # Return empty plot
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No significant genes found") +
             theme_void())
  }

  # Select genes based on user preference
  if (show_all) {
    display_genes <- sig_genes
    message("Showing all ", nrow(display_genes), " significant genes")
  } else {
    # Show top N genes from each group
    top_up <- sig_genes %>%
      dplyr::filter(!!as.name(significance_col) == "Upregulated") %>%
      dplyr::arrange(!!as.name(padj_col), dplyr::desc(abs(!!as.name(log2fc_col)))) %>%
      dplyr::slice_head(n = min(top_n_up, sum(sig_genes[[significance_col]] == "Upregulated")))

    top_down <- sig_genes %>%
      dplyr::filter(!!as.name(significance_col) == "Downregulated") %>%
      dplyr::arrange(!!as.name(padj_col), dplyr::desc(abs(!!as.name(log2fc_col)))) %>%
      dplyr::slice_head(n = min(top_n_down, sum(sig_genes[[significance_col]] == "Downregulated")))

    display_genes <- dplyr::bind_rows(top_up, top_down)
    message("Showing ", nrow(display_genes), " genes (",
            nrow(top_up), " up, ", nrow(top_down), " down)")
  }

  # Calculate absolute fold change for radius
  abs_fc <- abs(display_genes[[log2fc_col]])

  # Calculate significance for point size
  sig_values <- -log10(display_genes[[padj_col]])
  sig_values[is.na(sig_values)] <- 0
  sig_values[is.infinite(sig_values)] <- max(sig_values[is.finite(sig_values)], na.rm = TRUE) * 1.1

  # Prepare data for plotting
  plot_data <- display_genes %>%
    dplyr::mutate(
      # Assign to yin (down) or yang (up)
      yin_yang = ifelse(!!as.name(significance_col) == "Downregulated", "yin", "yang"),

      # Calculate angle for placement
      angle = ifelse(
        yin_yang == "yin",
        # Left half: 100° to 260° (5π/9 to 13π/9)
        (5*pi/9) + ((dplyr::row_number() - 1) / max(1, sum(yin_yang == "yin") - 1)) * (8*pi/9),
        # Right half: -80° to 80° (-4π/9 to 4π/9)
        (-4*pi/9) + ((dplyr::row_number() - 1) / max(1, sum(yin_yang == "yang") - 1)) * (8*pi/9)
      ),

      # Add slight random variation
      angle = angle + runif(n(), -0.05, 0.05),

      # Ensure genes stay in correct half
      angle = ifelse(
        yin_yang == "yin",
        pmin(pmax(angle, pi/2 + 0.1), 3*pi/2 - 0.1),
        pmin(pmax(angle, -pi/2 + 0.1), pi/2 - 0.1)
      ),

      # Radius based on absolute fold change
      raw_radius = scales::rescale(
        if (radius_power != 1) abs(!!as.name(log2fc_col))^radius_power else abs(!!as.name(log2fc_col)),
        to = c(min_radius, max_radius),
        from = if (radius_power != 1) {
          range(ifelse(is.na(abs(!!as.name(log2fc_col))), 0,
                       abs(!!as.name(log2fc_col))^radius_power), na.rm = TRUE)
        } else {
          range(abs(!!as.name(log2fc_col)), na.rm = TRUE)
        }
      ),

      radius = pmax(min_radius, pmin(max_radius, raw_radius)),

      # Calculate coordinates
      x = radius * cos(angle),
      y = radius * sin(angle),

      # Point size based on significance
      point_size = scales::rescale(
        sig_values,
        to = point_size_range,
        from = range(sig_values, na.rm = TRUE)
      ),

      # Gene label
      gene_label = dplyr::case_when(
        !is.na(!!as.name(gene_symbol_col)) &
          !!as.name(gene_symbol_col) != !!as.name(gene_id_col) &
          !!as.name(gene_symbol_col) != "" ~ as.character(!!as.name(gene_symbol_col)),
        TRUE ~ stringr::str_sub(!!as.name(gene_id_col), 1, 10)
      )
    )

  # Create background function
  create_background <- function(green_alpha = background_alpha, red_alpha = background_alpha) {
    # Generate points for yin semicircle (left side)
    yin_angle <- seq(pi/2, 3*pi/2, length.out = 100)
    yin_points <- data.frame(
      x = cos(yin_angle),
      y = sin(yin_angle)
    )

    # Generate points for yang semicircle (right side)
    yang_angle <- seq(-pi/2, pi/2, length.out = 100)
    yang_points <- data.frame(
      x = cos(yang_angle),
      y = sin(yang_angle)
    )

    # Create base plot
    base_plot <- ggplot() +

      # Outer circle
      annotate("path",
               x = cos(seq(0, 2*pi, length.out = 200)),
               y = sin(seq(0, 2*pi, length.out = 200)),
               color = "gray40", size = 1, fill = NA) +

      # Yin semicircle background
      annotate("polygon",
               x = c(yin_points$x, 0),
               y = c(yin_points$y, 0),
               fill = yin_color, alpha = green_alpha) +

      # Yang semicircle background
      annotate("polygon",
               x = c(yang_points$x, 0),
               y = c(yang_points$y, 0),
               fill = yang_color, alpha = red_alpha) +

      # Labels for yin and yang
      annotate("text",
               x = -0.7, y = 1.2,
               label = paste("YIN\n(Down:", sum(plot_data$yin_yang == "yin"), " genes)"),
               color = "darkgreen", size = 4, fontface = "bold", hjust = 0.5) +

      annotate("text",
               x = 0.7, y = 1.2,
               label = paste("YANG\n(Up:", sum(plot_data$yin_yang == "yang"), " genes)"),
               color = "darkred", size = 4, fontface = "bold", hjust = 0.5)

    return(base_plot)
  }

  # Generate subtitle if not provided
  if (is.null(subtitle)) {
    subtitle <- if (show_all) {
      paste("All significant genes (n =", nrow(plot_data), ")\n",
            "Radius = |log2FC|, Point size = -log10(padj)")
    } else {
      paste("Top", top_n_up, "up &", top_n_down, "down genes\n",
            "Radius = |log2FC|, Point size = -log10(padj)")
    }
  }

  # Create the plot
  p <- create_background() +

    # Add gene points
    geom_point(data = plot_data,
               aes(x = x, y = y,
                   size = point_size,
                   fill = yin_yang),
               shape = 21, color = "white", stroke = 0.7, alpha = 0.9) +

    # Color scale for points
    scale_fill_manual(
      values = c("yin" = point_yin_color, "yang" = point_yang_color),
      name = "Gene Regulation",
      breaks = c("yin", "yang"),
      labels = c("yin" = "Down-regulated", "yang" = "Up-regulated")
    ) +

    # Size scale
    scale_size_continuous(
      range = point_size_range,
      name = "Significance\n-log10(padj)",
      breaks = scales::rescale(
        quantile(sig_values, probs = c(0.1, 0.5, 0.9), na.rm = TRUE),
        to = point_size_range,
        from = range(sig_values, na.rm = TRUE)
      ),
      labels = c("Low", "Medium", "High")
    ) +

    # Title and subtitle
    labs(
      title = title,
      subtitle = subtitle,
      caption = if (!show_all && nrow(display_genes) < nrow(sig_genes)) {
        paste("Selected from", nrow(sig_genes), "total significant genes")
      } else {
        "Distance from center = |log2FoldChange|"
      }
    )

  # Label top N genes by significance
  top_n_to_label <- min(label_top_n, nrow(plot_data))
  top_labels <- plot_data %>%
    dplyr::arrange(!!as.name(padj_col), dplyr::desc(abs(!!as.name(log2fc_col)))) %>%
    dplyr::slice_head(n = top_n_to_label)

  message("Labeling top ", nrow(top_labels), " genes by significance")

  # Auto-calculate label size if not specified
  if (is.null(gene_label_size)) {
    # Try to get current device dimensions
    current_dev <- grDevices::dev.list()
    if (length(current_dev) > 0) {
      tryCatch({
        dev_size <- grDevices::dev.size()
        plot_device_width <- dev_size[1]
        plot_device_height <- dev_size[2]
      }, error = function(e) {
        plot_device_width <- 7
        plot_device_height <- 7
      })
    } else if (is.null(plot_device_width) || is.null(plot_device_height)) {
      plot_device_width <- 7
      plot_device_height <- 7
    }

    # Base size calculation
    base_size <- 4.5
    n_labels <- nrow(top_labels)

    # Adjust based on number of labels
    if (n_labels > 50) {
      gene_label_size <- base_size * 0.5
    } else if (n_labels > 30) {
      gene_label_size <- base_size * 0.7
    } else if (n_labels > 20) {
      gene_label_size <- base_size * 0.8
    } else if (n_labels > 10) {
      gene_label_size <- base_size * 0.9
    } else {
      gene_label_size <- base_size
    }

    # Adjust based on plot density
    plot_area <- plot_device_width * plot_device_height
    if (plot_area > 0) {
      density_factor <- nrow(plot_data) / plot_area
      gene_label_size <- gene_label_size / sqrt(1 + density_factor)
    }

    # Ensure reasonable bounds
    gene_label_size <- max(2, min(6, gene_label_size))
    message("Auto-calculated label size: ", round(gene_label_size, 2))
  }

  # Add labels
  p <- p +
    ggrepel::geom_label_repel(
      data = top_labels,
      aes(x = x, y = y, label = gene_label),
      size = gene_label_size,
      fill = ifelse(top_labels$yin_yang == "yin",
                    scales::alpha(yin_color, 0.85),
                    scales::alpha(yang_color, 0.85)),
      color = "black",
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.35,
      segment.size = 0.3,
      segment.color = "grey40",
      min.segment.length = 0.2,
      max.overlaps = Inf,
      max.time = 2,
      max.iter = 100000
    )

  # Add concentric circles to show radius scale
  p <- p +
    # Inner circle
    annotate("path",
             x = min_radius * cos(seq(0, 2*pi, length.out = 100)),
             y = min_radius * sin(seq(0, 2*pi, length.out = 100)),
             color = "gray70", size = 0.3, linetype = "dotted", alpha = 0.5) +

    # Middle circle
    annotate("path",
             x = (min_radius + max_radius)/2 * cos(seq(0, 2*pi, length.out = 100)),
             y = (min_radius + max_radius)/2 * sin(seq(0, 2*pi, length.out = 100)),
             color = "gray60", size = 0.3, linetype = "dashed", alpha = 0.5) +

    # Outer circle
    annotate("path",
             x = max_radius * cos(seq(0, 2*pi, length.out = 100)),
             y = max_radius * sin(seq(0, 2*pi, length.out = 100)),
             color = "gray40", size = 0.5, alpha = 0.7)

  # Theme customization
  p <- p +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold",
                                margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 11,
                                   margin = margin(b = 15)),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    coord_equal() +
    xlim(-1.2, 1.2) +
    ylim(-1.2, 1.2)

  # Add border if requested
  if (!is.null(border_color)) {
    p <- p +
      theme(panel.border = element_rect(color = border_color,
                                        fill = NA,
                                        size = border_size))
  }

  # Add legend
  if (show_legend) {
    p <- p +
      theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(20, 150, 20, 20)
      ) +
      guides(
        fill = guide_legend(
          order = 1,
          title = "Gene Regulation",
          override.aes = list(
            size = 5,
            alpha = 1,
            shape = 21,
            color = "white",
            stroke = 0.7
          )
        ),
        size = guide_legend(
          order = 2,
          title = "Significance\n-log10(padj)",
          override.aes = list(
            fill = "gray70",
            shape = 21,
            color = "white",
            stroke = 0.7,
            alpha = 0.9
          )
        )
      )
  } else {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}














#' Create Clean Pirate Plots for Gene Expression Visualization
#'
#' @description
#' Generate clean, minimalistic pirate plots for visualizing expression patterns
#' without statistical annotations. Focuses on clean visualization of expression
#' distributions across conditions.
#'
#' @import ggplot2
#' @import dplyr
#' @import ggbeeswarm
#' @import patchwork
#' @import DESeq2
#' @importFrom DESeq2 colData assay counts rlog varianceStabilizingTransformation
#' @importFrom utils head
#' @export
yy_pirateplot <- function(
    dds_object,
    results_df,
    gene_list,
    condition_col = "condition",
    group_names = NULL,
    group_colors = NULL,
    gene_names = NULL,
    plot_style = "pirate",
    normalization_method = "rlog",
    y_axis_title = "Normalized expression",
    title_size = 16,
    axis_title_size = 14,
    axis_text_size = 12,
    point_size = 3.5,
    violin_alpha = 0.6,
    beeswarm_cex = 3.5,
    layout_strategy = "auto",
    max_columns = 5,
    max_rows = 8,
    arrange_by = "input",
    save_plots = FALSE,
    output_dir = "YinYangPlot_PiratePlots_Clean",
    plot_width = 16,
    plot_height = 12,
    dpi = 300,
    return_data = TRUE,
    seed = 42
) {

  set.seed(seed)

  # ============================================================================
  # 1. VALIDATE INPUTS AND SETUP
  # ============================================================================

  required_packages <- c("ggplot2", "dplyr", "ggbeeswarm", "patchwork", "DESeq2")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Required package(s) not installed: ", paste(missing_packages, collapse = ", "))
  }

  if (!inherits(dds_object, "DESeqDataSet")) {
    stop("'dds_object' must be a DESeqDataSet object")
  }

  if (!is.data.frame(results_df)) {
    stop("'results_df' must be a data frame")
  }

  if (!is.character(gene_list) || length(gene_list) == 0) {
    stop("'gene_list' must be a non-empty character vector")
  }

  col_data <- colData(dds_object)

  if (!condition_col %in% colnames(col_data)) {
    stop("Column '", condition_col, "' not found in colData(dds_object)")
  }

  original_groups <- unique(as.character(col_data[[condition_col]]))
  n_groups <- length(original_groups)

  if (n_groups < 2) {
    stop("Need at least 2 groups for comparison.")
  }

  if (is.null(group_names)) {
    group_names <- stats::setNames(original_groups, original_groups)
  }

  display_groups <- group_names[original_groups]
  names(display_groups) <- original_groups

  if (is.null(group_colors)) {
    palette <- c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#999999"
    )[seq_len(n_groups)]
    group_colors <- stats::setNames(palette, original_groups)
  }

  display_colors <- group_colors[original_groups]
  names(display_colors) <- display_groups

  # ============================================================================
  # 2. PREPARE EXPRESSION DATA
  # ============================================================================

  if (normalization_method == "rlog") {
    norm_data <- DESeq2::rlog(dds_object, blind = FALSE)
    norm_values <- assay(norm_data)
  } else if (normalization_method == "vst") {
    norm_data <- DESeq2::varianceStabilizingTransformation(dds_object, blind = FALSE)
    norm_values <- assay(norm_data)
  } else {
    norm_counts <- DESeq2::counts(dds_object, normalized = TRUE)
    norm_values <- log2(norm_counts + 1)
  }

  # ============================================================================
  # 3. MATCH GENES AND PREPARE PLOT DATA
  # ============================================================================

  plot_data_all <- data.frame()
  gene_info <- data.frame()
  not_found <- character()

  for (i in seq_along(gene_list)) {
    gene_query <- gene_list[i]
    gene_display <- ifelse(!is.null(gene_names) && length(gene_names) >= i,
                           gene_names[i], gene_query)

    found <- FALSE

    if ("gene_symbol" %in% colnames(results_df)) {
      matches <- results_df[results_df$gene_symbol == gene_query, ]
      if (nrow(matches) > 0) {
        gene_row <- matches[1, ]
        found <- TRUE
      }
    }

    if (!found && "gene_id" %in% colnames(results_df)) {
      matches <- results_df[grep(gene_query, results_df$gene_id), ]
      if (nrow(matches) > 0) {
        gene_row <- matches[1, ]
        found <- TRUE
      }
    }

    if (!found) {
      not_found <- c(not_found, gene_query)
      next
    }

    gene_id_found <- gene_row$gene_id

    if (!gene_id_found %in% rownames(norm_values)) {
      not_found <- c(not_found, gene_query)
      next
    }

    expression_values <- norm_values[gene_id_found, ]

    gene_df <- data.frame(
      Gene = gene_display,
      Sample = names(expression_values),
      Expression = as.numeric(expression_values),
      Display_Group = display_groups[
        as.character(col_data[names(expression_values), condition_col])
      ],
      stringsAsFactors = FALSE
    )

    plot_data_all <- rbind(plot_data_all, gene_df)
  }

  if (nrow(plot_data_all) == 0) {
    stop("No genes found.")
  }

  plot_data_all$Display_Group <- factor(
    plot_data_all$Display_Group,
    levels = display_groups
  )

  # ============================================================================
  # 5. CREATE CLEAN INDIVIDUAL PLOTS
  # ============================================================================

  individual_plots <- list()
  unique_genes <- unique(plot_data_all$Gene)

  for (gene in unique_genes) {

    gene_data <- plot_data_all[plot_data_all$Gene == gene, ]

    y_min <- min(gene_data$Expression) * 0.98
    y_max <- max(gene_data$Expression) * 1.02

    p <- ggplot(gene_data, aes(x = Display_Group, y = Expression, fill = Display_Group)) +
      geom_violin(trim = TRUE, alpha = violin_alpha, width = 0.7) +
      geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm(
        size = point_size,
        alpha = 0.8,
        shape = 21,
        color = "black",
        cex = beeswarm_cex
      ) +
      scale_fill_manual(values = display_colors) +
      labs(title = gene, x = "", y = y_axis_title) +
      theme_minimal(base_size = axis_text_size) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = title_size
        ),
        axis.text.x = element_text(
          size = axis_text_size,
          face = "bold",
          color = "black"
        ),
        axis.text.y = element_text(size = axis_text_size),
        legend.position = "none",

        # REMOVE GRID LINES
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
      ) +
      ylim(y_min, y_max)

    individual_plots[[gene]] <- p
  }

  combined_plot <- patchwork::wrap_plots(individual_plots)

  result <- list(
    individual_plots = individual_plots,
    combined_plot = combined_plot,
    plot_data = if (return_data) plot_data_all else NULL,
    not_found = not_found
  )

  invisible(result)
}

# ==============================================================================
# SIMPLIFIED VERSION
# ==============================================================================

#' Quick Clean Pirate Plot
#' @export
yy_pirateplot_simple <- function(dds_object, results_df, gene, ...) {
  result <- yy_pirateplot(
    dds_object = dds_object,
    results_df = results_df,
    gene_list = gene,
    save_plots = FALSE,
    return_data = FALSE,
    ...
  )
  result$individual_plots[[1]]
}



















