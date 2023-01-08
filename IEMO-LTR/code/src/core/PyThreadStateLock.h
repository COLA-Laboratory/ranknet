//
// Created by gylai on 2021/8/2.
//

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>

class PyThreadStateLock
        {
        public:
            PyThreadStateLock()
            {
                state = PyGILState_Ensure();//get GIL
            }

            ~PyThreadStateLock()
            {
                PyGILState_Release( state );
            }
        private:
            PyGILState_STATE state;
        };