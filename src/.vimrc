let g:ale_linters = {'cpp': ['clang']}  " use clang compiler
let cpp_flags = '-std=c++14 -Wall '  " c++ standard & show all warnings
let cpp_flags = cpp_flags . '-I.' . ' '  " for headres
let cpp_flags = cpp_flags . '-I../extern/eigen3' . ' '  " for eigen
let cpp_flags = cpp_flags . '-I../extern/pybind11/include' . ' '  " for eigen
let cpp_flags = cpp_flags . '-I/usr/local/include/python3.7m' . ' '  " for eigen
let cpp_flags = cpp_flags . '-Wno-unknown-warning-option' . ' '  " supress eigen warnings
let g:ale_cpp_cc_options = cpp_flags  " set the flag for linting
let g:ale_cpp_clang_options = cpp_flags  " set the flag for linting


" use ctrl+j/k to jump between errors
nmap <silent> <C-k> <Plug>(ale_previous_wrap)
nmap <silent> <C-j> <Plug>(ale_next_wrap)
