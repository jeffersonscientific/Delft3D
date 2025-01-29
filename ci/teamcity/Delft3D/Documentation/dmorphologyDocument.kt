package Delft3D.windows

import jetbrains.buildServer.configs.kotlin.*
import jetbrains.buildServer.configs.kotlin.buildFeatures.*
import jetbrains.buildServer.configs.kotlin.buildSteps.*
import jetbrains.buildServer.configs.kotlin.failureConditions.*

import Delft3D.template.*

object DMorphologyValidationDocument : BuildType({
    templates(TemplateDocumentationBuild)

    name = "D-Morphology - Validation and Functionality document (Latex/PDF)"

    params {
        param("engine_dir", "e02_dflowfm")
        param("engine_name", "dmorphology")
    }
})
