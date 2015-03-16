(TeX-add-style-hook
 "p2"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("subfiles" "../main.tex")))
   (TeX-run-style-hooks
    "latex2e"
    "subfiles"
    "subfiles10")
   (LaTeX-add-labels
    "cha:problem-2"
    "sec:introduction"
    "eq:w"
    "eq:P"
    "eq:C"
    "eq:planeStress"
    "sec:consistency"
    "sec:mfi"
    "sec:isotropy"
    "sec:NR"
    "sec:calc"
    "fig:neoHcon"
    "fig:psneoHcon"
    "fig:uniaxial"
    "fig:equibiaxial")))

