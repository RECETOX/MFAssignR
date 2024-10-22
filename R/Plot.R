create_mass_spectrum_plot <- function(records1, unassigned, rawpeaks) {
    plot <- ggplot2::ggplot() + ggplot2::geom_segment(data = records1, size = 0.7, ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "RA"), color = "green") +
        ggplot2::geom_segment(data = records1, size = 0.7, ggplot2::aes_string(x = "C13_mass", xend = "C13_mass", y = 0, yend = "C13_Abund"), color = "blue") + ggplot2::geom_segment(data = records1,
        size = 0.7, ggplot2::aes_string(x = "C13_mass2", xend = "C13_mass2", y = 0, yend = "C13_Abund2"), color = "blue") + ggplot2::geom_segment(data = records1, size = 0.7, ggplot2::aes_string(x = "S34_mass",
        xend = "S34_mass", y = 0, yend = "S34_Abund"), color = "blue") + ggplot2::geom_segment(data = unassigned, size = 0.7, ggplot2::aes_string(x = "mass", xend = "mass", y = 0,
        yend = "RA"), color = "red") + ggplot2::coord_cartesian(xlim = c(min(rawpeaks$mass), max(rawpeaks$mass))) + ggplot2::theme_bw() + ggplot2::labs(x = "Ion Mass", y = "Abundance",
        title = "Assignment Mass Spectrum", color = "DBE") + ggplot2::theme(axis.title = ggplot2::element_text(size = 15, face = "bold"), strip.text = ggplot2::element_text(size = 15,
        face = "bold"), axis.text = ggplot2::element_text(size = 15, face = "bold"), legend.title = ggplot2::element_text(face = "bold", size = 15), legend.text = ggplot2::element_text(face = "bold",
        size = 15), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16,
        face = "bold"))
    return(plot)
}

create_error_plot <- function(Unambig, Ambigout, records1) {
    plot <- ggplot2::ggplot() + ggplot2::geom_point(data = Unambig, ggplot2::aes_string(x = "Exp_mass", y = "AE_ppm", color = "Tag"), alpha = 1/3) + ggplot2::geom_point(data = Ambigout,
        ggplot2::aes_string(x = "Exp_mass", y = "AE_ppm", color = "Tag"), alpha = 1/3) + ggplot2::coord_cartesian(xlim = c(min(records1$Exp_mass), max(records1$Exp_mass)), ylim = c(min(records1$AE_ppm),
        max(records1$AE_ppm))) + ggplot2::scale_colour_manual(name = "Ambiguity", values = c(Unambiguous = "blue", Ambiguous = "red")) + ggplot2::labs(x = "Ion Mass", y = "Absolute Error (ppm)",
        color = "Ambiguity", title = "Error Plot") + ggplot2::theme_bw() + ggplot2::theme(axis.title = ggplot2::element_text(size = 15, face = "bold"), strip.text = ggplot2::element_text(size = 15,
        face = "bold"), axis.text = ggplot2::element_text(size = 15, face = "bold"), legend.title = ggplot2::element_text(face = "bold", size = 11), legend.text = ggplot2::element_text(face = "bold",
        size = 10), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16,
        face = "bold"))
    return(plot)
}

create_mass_spectrum_group_plot <- function(PD, form_palette) {
    plot <- ggplot2::ggplot() + ggplot2::geom_segment(data = PD, size = 0.7, ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "RA", color = "Tag2")) + ggplot2::facet_wrap(~Tag,
        ncol = 1, scales = "free_y") + ggplot2::scale_colour_manual(name = "Groups", values = form_palette) + ggplot2::theme_bw() + ggplot2::labs(x = "Ion Mass", y = "Abundance",
        title = "Assignment Mass Spectrum", color = "DBE") + ggplot2::theme(axis.title = ggplot2::element_text(size = 15, face = "bold"), strip.text = ggplot2::element_text(size = 15,
        face = "bold"), axis.text = ggplot2::element_text(size = 15, face = "bold"), legend.title = ggplot2::element_text(face = "bold", size = 12), legend.text = ggplot2::element_text(face = "bold",
        size = 12), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16,
        face = "bold"))
    return(plot)
}

create_van_krevelen_plot <- function(PD, form_palette) {
    plot <- ggplot2::ggplot() + ggplot2::geom_point(data = PD, ggplot2::aes_string(x = "O_C", y = "H_C", color = "Tag2"), alpha = 1/3) + ggplot2::facet_wrap(~Tag, ncol = 2) + ggplot2::scale_colour_manual(name = "Groups",
        values = form_palette) + ggplot2::labs(x = "Oxygen-to-Carbon Ratio", y = "Hydrogen-to-Carbon Ratio", color = "Groups", title = "van Krevelen Plot") + ggplot2::theme_bw() +
        ggplot2::theme(axis.title = ggplot2::element_text(size = 15, face = "bold"), strip.text = ggplot2::element_text(size = 15, face = "bold"), axis.text = ggplot2::element_text(size = 15,
            face = "bold"), legend.title = ggplot2::element_text(face = "bold", size = 11), legend.text = ggplot2::element_text(face = "bold", size = 10), panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16, face = "bold"))
    return(plot)
}
