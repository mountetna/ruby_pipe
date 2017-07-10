require_relative 'example/script'

describe "Pipeline::Script" do
  it "returns a list of tasks" do
    script = MutationDetection.new

    expect(script.tasks).to be([])
  end
end
