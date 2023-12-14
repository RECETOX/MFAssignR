#' Assigns all possible combinations of C, H, and O for a given MW
#'
#' FindCoreFormulae assigns all possible combinations of C, H, and O for a
#' given MW and subject to constraints in the \code{\link{ValidFormula}} function.
#' This is a subfunction of the MFAssign function and will not work
#' independently of that function.
#'
#' @param env environment from \code{\link{MFAssign}} function
#'
#' @examples
#' env <- environment()
#' FindCoreFormulae(env)
#'
#' @export
FindCoreFormulae2_Halo <- function(env) {
  if (Even(env$coreRNM) == FALSE) {
    return()
  }

  Output <- create_output(env)

  env$records <- list(
    RA = (env$RA), coreNM = env$coreRNM, Exp_mass = env$ionEM,
    C = Output$C,
    H = Output$H, O = Output$O,
    N = env$loop[CompFactorToInt2("N")], S = env$loop[CompFactorToInt2("S")],
    P = env$loop[CompFactorToInt2("P")], E = env$loop[CompFactorToInt2("E")],
    S34 = env$loop[CompFactorToInt2("S34")], N15 = env$loop[CompFactorToInt2("N15")],
    D = env$loop[CompFactorToInt2("D")], Cl = env$loop[CompFactorToInt2("Cl")],
    Fl = env$loop[CompFactorToInt2("Fl")],
    Cl37 = env$loop[CompFactorToInt2("Cl37")], Br = env$loop[CompFactorToInt2("Br")],
    Br81 = env$loop[CompFactorToInt2("Br81")], I = env$loop[CompFactorToInt2("I")],
    M = env$loop[CompFactorToInt2("M")],
    NH4 = env$loop[CompFactorToInt2("NH4")], POE = env$loop[CompFactorToInt2("POE")],
    NOE = env$loop[CompFactorToInt2("NOE")],
    Z = env$loop[CompFactorToInt2("Z")],
    Neutral_mass = Output$XEM, CHO_mass = Output$CEM,
    CHO_Err = Output$xemErr, Ratio = Output$Ratio2
  )
}
