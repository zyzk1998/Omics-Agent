!macro customInstall
  nsExec::ExecToStack 'taskkill /F /IM "Omics Agent.exe" /T'
  nsExec::ExecToStack 'taskkill /F /IM local_sidecar.exe /T'
  Sleep 1000
!macroend
