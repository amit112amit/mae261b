(TeX-add-style-hook
 "p1"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("subfiles" "../main.tex")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "subfiles"
    "subfiles10")
   (LaTeX-add-labels
    "fig:linDom"
    "eq:bary"
    "fig:quadDom"
    "tab:gaussQuad"
    "fig:linCon"
    "fig:quadCon"
    "fig:consistency"
    "tab:rank")))

