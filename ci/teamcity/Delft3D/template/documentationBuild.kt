package Delft3D.template

import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildSteps.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.*
import jetbrains.buildServer.configs.kotlin.triggers.*


object TemplateDocumentationBuild : Template({
    name = "Generate validation and functionality reports"
    description = "This build configuration generates validation and functionality reports for the Delft3D engine."
    buildNumberPattern = "%build.vcs.number%"

    artifactRules = """
        %engine_dir%/*.log=>logging
        %engine_dir%/doc/validation/*.pdf=>pdf
        %engine_dir%/doc/validation/*.log=>logging
        %engine_dir%/doc/functionalities/*.pdf=>pdf
        %engine_dir%/doc/functionalities/*.log=>logging
        %engine_dir%/*/doc/*.pdf=>pdf/functionality
        %engine_dir%/*/doc/*.log=>logging/functionality
    """.trimIndent()

    params {
        param("s3_dsctestbench_accesskey", DslContext.getParameter("s3_dsctestbench_accesskey"))
        password("s3_dsctestbench_secret", "credentialsJSON:7e8a3aa7-76e9-4211-a72e-a3825ad1a160")
    }

    vcs {
        root(DslContext.settingsRoot)
        cleanCheckout = true
    }

    steps {
        python {
            name = "Checkout Testbench cases from MinIO"
            id = "CHECKOUT_TESTBENCH_CASES_FROM_MINIO"
            environment = venv {
            requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            command = file {
                filename = "test/deltares_testbench/scripts/download_from_s3.py"
                scriptArguments = "--access_key %s3_dsctestbench_accesskey% --secret_key %s3_dsctestbench_secret% --engine_dir %engine_dir%"
            }
            }
        }
        python {
            name = "Checkout Testbench cases from MinIO"
            id = "CHECKOUT_TESTBENCH_CASES_FROM_MINIO"
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            command = file {
                filename = "test/deltares_testbench/scripts/checkout_testbench_cases.py"
            }
        }
        python {
            name = "Update infrastructure for functionality report"
            id = "UPDATE_INFRASTRUCTURE_FUNCTIONALITY_REPORT"
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            command = file {
                filename = "test/deltares_testbench/scripts/update_functionality_report.py"
                scriptArguments = "--reldir ./%engine_dir%"
            }
        }
        python {
            name = "Generate report"
            id = "GENERATE_REPORT"
            environment = venv {
                requirementsFile = "test/deltares_testbench/scripts/requirements.txt"
            }
            command = file {
                filename = "test/deltares_testbench/scripts/generate_report.py"
                scriptArguments = "--texfile %engine_dir%/doc/validation/%engine_name%_validation_doc.tex"
            }
        }
    }

    requirements {
        startsWith("teamcity.agent.jvm.os.name", "Windows", "RQ_2470")
    }

    triggers {
        finishBuildTrigger {
            buildType = "Dimr_DimrCollector"
            branchFilter = "+:<default>"
        }
    }
})
