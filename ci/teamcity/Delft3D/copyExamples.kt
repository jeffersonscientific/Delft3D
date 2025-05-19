import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.*
import jetbrains.buildServer.configs.kotlin.buildSteps.*
import jetbrains.buildServer.configs.kotlin.triggers.*
import Delft3D.template.*
import Delft3D.step.*
import Delft3D.linux.*
import Delft3D.windows.*

object CopyExamples : BuildType({
    id("DHydro_ExampleCases_CopyDelft3dfmExampleCases")
    name = "Copy delft3dfm example cases"
    description = "Copy example files to P drive"
    buildNumberPattern = "%build.vcs.number%"

    templates(
        TemplateMonitorPerformance
    )

    vcs {
        root(DslContext.settingsRoot)
    }

    params {
        param("DEST_DIR", """\\directory.intra\PROJECT\d-hydro\dimrset\examples""")
    }

    steps {
        python {
            name = "Run TestBench.py"
            workingDir = "test/deltares_testbench/"
            environment = venv {
                requirementsFile = "pip/win-requirements.txt"
            }
            command = file {
                filename = "TestBench.py"
                scriptArguments = """
                    --username "%s3_dsctestbench_accesskey%"
                    --password "%s3_dsctestbench_secret%"
                    --reference
                    --config "configs/singularity/dimr/dimr_smoke_test_lnx64.xml"
                    --filter "testcase=e02_f014_c001,e109_f01_c010,e112_f01_c11"
                    --skip-run
                    --skip-post-processing
                    --log-level DEBUG
                    --parallel
                    --teamcity
                """.trimIndent()
            }
        }
        python {
            name = "Copy example files to P drive"
            environment = venv {
                requirementsFile = ""
                pipArgs = "--editable ./ci/python"
            }
            command = module {
                module = "ci_tools.example_utils.copy_examples"
                scriptArguments = """
                    --dest_dir
                    %DEST_DIR% 
                    --tc_logging
                """.trimIndent()
            }
        }
    }

    if (DslContext.getParameter("environment") == "production") {
        triggers {
            vcs {
                branchFilter = "+:<default>"
            }
        }
    }
    requirements {
        contains("teamcity.agent.jvm.os.name", "Windows")
    }
})