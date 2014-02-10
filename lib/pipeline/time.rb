module TimeFormat
  def human_time secs
    case secs
    when 0..59
      "#{secs} seconds"
    when 60..599
      "#{secs/60} minutes and #{secs % 60} seconds"
    when 599..3599
      "#{secs/60} minutes"
    when 3600..35999
      "#{secs/3600} hours and #{(secs%3600)/60} minutes"
    else
      "#{secs/3600} hours"
    end
  end
end
