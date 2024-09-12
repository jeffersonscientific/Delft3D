import jetbrains.buildServer.configs.kotlin.*

import Testbench.*

version = "2024.03"

project {
    description = "contact: BlackOps (black-ops@deltares.nl)"

    buildType(LinuxAll)

    buildTypesOrder = arrayListOf(LinuxAll)
}
