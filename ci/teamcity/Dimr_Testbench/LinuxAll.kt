package testbench

import jetbrains.buildServer.configs.kotlin.*

object LinuxAll : BuildType({
    templates(LinuxTestTemplate)
    name = "Linux (all)"

    params {
        param("branch", "all")
    }
})