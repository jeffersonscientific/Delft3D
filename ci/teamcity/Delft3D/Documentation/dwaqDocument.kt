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
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            name = "Generate functionality report"
            id = "GENERATE_FUNCTIONALITY_REPORT"
            command = file {
                filename = "test/deltares_testbench/scripts/generate_functionality_report.py"
                scriptArguments = "--engine_dir_name %engine_dir%"
            }
        }
        stepsOrder = arrayListOf("CHECKOUT_TESTBENCH_CASES_FROM_MINIO", "UPDATE_INFRASTRUCTURE_FUNCTIONALITY_REPORT", "GENERATE_FUNCTIONALITY_REPORT", "GENERATE_REPORT")
    }
})
