requirejs.config({
    paths: {
        'mathjax': 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS_HTML',
    },
    shim: {
        'mathjax' : {
            exports: "MathJax"
        },
    }
});

require(['mathjax'], function(MathJax) {
    MathJax.Hub.Config({
        TeX: {
            Macros: {
                defd: "≝",
                abs: ["\\left|#1\\right|",1],
                vec: ["\\mathbf{#1}",1],
                mat: ["\\mathsf{#1}",1],
                conj: ["#1^*",1],
                ce: "\\mathrm{e}",
                im: "\\mathrm{i}",
                diff : ["\\mathrm{d}#1\\,",1],
                Beta : "\\mathrm{B}", // You've gotta be kidding
                divdiff : "⍋",
                bmat: ["\\begin{bmatrix}#1\\end{bmatrix}",1]
            }
        }
    });
})


