package Delft3D.windows

import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.*
import jetbrains.buildServer.configs.kotlin.buildSteps.*
import jetbrains.buildServer.configs.kotlin.failureConditions.*


import Delft3D.template.*

object DflowfmDwavesValidationDocument : BuildType({
    templates(TemplateDocumentationBuild)

    name = "D-Flow FM, D-Waves - Validation and Functionality document (Latex/PDF)"

    steps {
        python {
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            name = "Generate functionality report"
            command = file {
                filename = "test/deltares_testbench/scripts/generate_functionality_report.py"
                scriptArguments = "--engine_dir_name %engine_dir%"
            }
        }
    }

    params {
        param("engine_dir", "e100_dflowfm-dwaves")
        param("engine_name", "dflowfm-dwaves")
    }
})
