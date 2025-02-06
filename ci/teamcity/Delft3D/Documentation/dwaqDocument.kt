package Delft3D.windows

import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.*
import jetbrains.buildServer.configs.kotlin.buildSteps.*
import jetbrains.buildServer.configs.kotlin.failureConditions.*

import Delft3D.template.*

object DWaqValidationDocument : BuildType({
    templates(TemplateDocumentationBuild)

    name = "D-Water Quality- Validation and Functionality document (Latex/PDF)"

    params {
        param("engine_dir", "e03_waq")
        param("engine_name", "dwaq")
    }

    steps {
        python {
            name = "Generate functionality report"
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            command = file {
                filename = "test/deltares_testbench/scripts/generate_functionality_report.py"
                scriptArguments = "--engine_dir_name %engine_dir%"
            }
        }
    }
})
